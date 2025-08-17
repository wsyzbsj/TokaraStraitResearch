import toml
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime, timedelta
import cartopy.io.shapereader as shpreader
from shapely.geometry import box
from shapely.strtree import STRtree
import matplotlib.path as mpath

# 设置
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
config = toml.load('config/config.toml')
# 局部放大区域
lon_min = config['zone_config']['lon_min_zoom']
lon_max = config['zone_config']['lon_max_zoom']
lat_min = config['zone_config']['lat_min_zoom']
lat_max = config['zone_config']['lat_max_zoom']

# 地球物理参数
OMEGA = 7.292115e-5  # 地球自转角速度 (rad/s)
R = 6371000.0        # 地球半径 (m)

# 读取流场数据
dataset = netCDF4.Dataset(r'D:\Documents\Data\HYCOM\uv3z_2024071709.nc4', 'r')
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
u = dataset.variables['water_u'][:]     # (time,depth,lat,lon)
v = dataset.variables['water_v'][:]     # (time,depth,lat,lon)
depth_levels = dataset.variables['depth'][:]  # 深度层数据
time_val = dataset.variables['time'][:][0]

# 筛选数据 - 确保使用相同的索引
lat_idx = np.where((lat >= lat_min) & (lat <= lat_max))[0]
lon_idx = np.where((lon >= lon_min) & (lon <= lon_max))[0]
lat_used = lat[lat_idx]
lon_used = lon[lon_idx]

# 获取深度层数量
n_depths = len(depth_levels)

# 创建规则网格
print('创建网格')
grid_x, grid_y = np.mgrid[lon_min:lon_max:200j, lat_min:lat_max:200j]

print('读取海岸线数据')
shp_path = r'D:\Documents\Data\海岸线数据\land_polygons.shp'
reader = shpreader.Reader(shp_path)
# 创建空间索引
geoms = list(reader.geometries())
tree = STRtree(geoms)
# 创建目标区域边界框
target_bbox = box(lon_min, lat_min, lon_max, lat_max)
# 查询相交的几何体
intersect_indices = tree.query(target_bbox, predicate="intersects")
valid_geoms = [geoms[i] for i in intersect_indices]

# 使用更高效的方法计算陆地掩膜
print('使用matplotlib.Path计算陆地掩膜')
land_mask = np.zeros(grid_x.shape, dtype=bool)
# 收集所有多边形的顶点
all_polygons = []
for geom in valid_geoms:
    if geom.geom_type == 'Polygon':
        all_polygons.append(geom)
    elif geom.geom_type == 'MultiPolygon':
        all_polygons.extend(geom.geoms)
# 创建包含所有多边形的复合路径
vertices = []
codes = []
for poly in all_polygons:
    # 处理外部环
    ext_vertices = np.array(poly.exterior.coords)
    vertices.append(ext_vertices)
    codes.append([mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(ext_vertices) - 2) + [mpath.Path.CLOSEPOLY])
    # 处理内部环（孔洞）
    for interior in poly.interiors:
        int_vertices = np.array(interior.coords)
        vertices.append(int_vertices)
        codes.append([mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(int_vertices) - 2) + [mpath.Path.CLOSEPOLY])
# 创建复合路径
compound_path = mpath.Path(np.concatenate(vertices), np.concatenate(codes))
# 一次性检查所有网格点
points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
inside = compound_path.contains_points(points)
land_mask = inside.reshape(grid_x.shape)

# 设置投影和范围
proj = ccrs.PlateCarree()
extents = [lon_min, lon_max, lat_min, lat_max]

time_str = (datetime(2000,1,1)+timedelta(hours=time_val)).strftime('%Y%m%d%H')

# 计算经度和纬度的网格间距（度）
dlon = np.mean(np.diff(lon_used))
dlat = np.mean(np.diff(lat_used))

# 计算层厚度（基于HYCOM深度层数据）
# 假设深度层表示层中心深度
# 计算层边界（层顶部和底部深度）
layer_boundaries = np.zeros(len(depth_levels) + 1)
# 第一层顶部为0（海面）
layer_boundaries[0] = 0
# 中间层边界为相邻层深度的平均值
for i in range(1, len(depth_levels)):
    layer_boundaries[i] = (depth_levels[i-1] + depth_levels[i]) / 2
# 最后一层底部为最后一层深度的1.5倍（或根据实际情况调整）
layer_boundaries[-1] = depth_levels[-1] * 1.5

# 计算每层厚度
layer_thicknesses = layer_boundaries[1:] - layer_boundaries[:-1]

# 创建原始数据网格点（用于插值）
lon_mesh, lat_mesh = np.meshgrid(lon_used, lat_used)
points_data = np.column_stack((lon_mesh.ravel(), lat_mesh.ravel()))

# 循环处理每个深度层
for depth_idx in range(n_depths):
    print(f"处理深度层 {depth_idx+1}/{n_depths} (深度: {depth_levels[depth_idx]} m)")

    # 获取当前层的流速数据
    u_layer = u[0, depth_idx, lat_idx[0]:lat_idx[-1]+1, lon_idx[0]:lon_idx[-1]+1]
    v_layer = v[0, depth_idx, lat_idx[0]:lat_idx[-1]+1, lon_idx[0]:lon_idx[-1]+1]

    # 计算科氏参数 f = 2Ωsinφ (二维数组)
    f = 2 * OMEGA * np.sin(np.deg2rad(lat_used))
    f_2d = np.tile(f[:, np.newaxis], (1, len(lon_used)))

    # 计算相对涡度 ζ = ∂v/∂x - ∂u/∂y
    # 初始化涡度场
    zeta = np.zeros_like(u_layer)

    # 计算网格间距（米）
    dx = R * np.cos(np.deg2rad(lat_used)) * np.deg2rad(dlon)  # 经度方向间距 (m)
    dy = R * np.deg2rad(dlat)  # 纬度方向间距 (m)

    # 计算经度方向的偏导数 (∂v/∂x)
    dv_dx = np.zeros_like(v_layer)
    for i in range(len(lat_used)):
        # 处理边界情况
        if i < len(lat_used) and i < dv_dx.shape[0]:
            if dv_dx.shape[1] > 2:  # 确保有足够的点进行中心差分
                dv_dx[i, 1:-1] = (v_layer[i, 2:] - v_layer[i, :-2]) / (2 * dx[i])
                # 边界处理
                dv_dx[i, 0] = (v_layer[i, 1] - v_layer[i, 0]) / dx[i]
                dv_dx[i, -1] = (v_layer[i, -1] - v_layer[i, -2]) / dx[i]
            else:  # 如果点数不足，使用前向/后向差分
                if dv_dx.shape[1] > 1:
                    dv_dx[i, 0] = (v_layer[i, 1] - v_layer[i, 0]) / dx[i]
                    dv_dx[i, -1] = (v_layer[i, -1] - v_layer[i, -2]) / dx[i]

    # 计算纬度方向的偏导数 (∂u/∂y)
    du_dy = np.zeros_like(u_layer)
    for j in range(len(lon_used)):
        # 处理边界情况
        if j < len(lon_used) and j < du_dy.shape[1]:
            if du_dy.shape[0] > 2:  # 确保有足够的点进行中心差分
                du_dy[1:-1, j] = (u_layer[2:, j] - u_layer[:-2, j]) / (2 * dy)
                # 边界处理
                du_dy[0, j] = (u_layer[1, j] - u_layer[0, j]) / dy
                du_dy[-1, j] = (u_layer[-1, j] - u_layer[-2, j]) / dy
            else:  # 如果点数不足，使用前向/后向差分
                if du_dy.shape[0] > 1:
                    du_dy[0, j] = (u_layer[1, j] - u_layer[0, j]) / dy
                    du_dy[-1, j] = (u_layer[-1, j] - u_layer[-2, j]) / dy

    # 计算相对涡度
    zeta = dv_dx - du_dy

    # 获取当前层厚度
    h_layer = np.ones_like(zeta) * layer_thicknesses[depth_idx]

    # 计算位涡 PV = (f + ζ) / h
    PV = (f_2d + zeta) / h_layer

    # 创建图形
    fig = plt.figure(figsize=(12, 10), dpi=150)
    ax = plt.axes(projection=proj)
    ax.set_extent(extents, crs=proj)
    ax.set_facecolor('lightgray')
    ax.set_title(f"{time_str} HYCOM 深度层 {depth_idx+1} ({depth_levels[depth_idx]} m) 位涡", fontsize=16)

    # 插值到位涡到规则网格
    grid_PV = griddata(points_data, PV.ravel(), (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 应用陆地掩膜
    grid_PV[land_mask] = np.nan

    # 创建自定义颜色映射
    ocean_cmap = plt.get_cmap('coolwarm').copy()

    # 计算位涡的对称范围
    pv_max = np.nanmax(np.abs(grid_PV))
    vmin, vmax = -pv_max, pv_max

    # 绘制位涡填色图
    cf = ax.contourf(grid_x, grid_y, grid_PV, levels=50, cmap=ocean_cmap,
                     transform=proj, vmin=vmin, vmax=vmax)

    # 添加colorbar
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
    cbar.set_label('位涡 (m⁻¹s⁻¹)', fontsize=12)

    # 绘制海岸线
    for geom in valid_geoms:
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            facecolor='none',
            edgecolor='black',
            zorder=5
        )

    # 添加网格线
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # 保存图像
    filename = f"output/HYCOM/位涡/PV_depth_{depth_idx+1}_{depth_levels[depth_idx]}m_{time_str}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"已保存: {filename}")

# 关闭数据集
dataset.close()
print("所有深度层处理完成")