import toml
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime, timedelta
import pandas as pd
import cartopy.io.shapereader as shpreader
from shapely.geometry import box
from shapely.prepared import prep
from shapely.strtree import STRtree
import matplotlib.path as mpath

# 设置
config = toml.load('config/config.toml')
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
# 局部放大区域
lon_min = config['zone_config']['lon_min_zoom']
lon_max = config['zone_config']['lon_max_zoom']
lat_min = config['zone_config']['lat_min_zoom']
lat_max = config['zone_config']['lat_max_zoom']

# 读取流场数据
dataset = netCDF4.Dataset(r'D:\Documents\code\Pycharm\FVCOM_analysis\kuroshio_2024_0717\kuroshio_2024_0717.nc', 'r')
nv = dataset.variables['nv'][:]
time = dataset.variables['time'][0]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
latc = dataset.variables['latc'][:]
lonc = dataset.variables['lonc'][:]
u = dataset.variables['u'][:]
v = dataset.variables['v'][:]
times = dataset.variables['Times'][:]
h = dataset.variables['h'][:]  # 水深
zeta = dataset.variables['zeta'][:]  # 水位
siglay = dataset.variables['siglay'][:]  # sigma层
temp = dataset.variables['temp'][:]

# 计算每个网格点的实际深度
print('计算每个网格点的深度')
depth_all = np.zeros((zeta.shape[0], siglay.shape[0], len(lat)), dtype=float)
for i_time in range(zeta.shape[0]):
    for i_lev in range(siglay.shape[0]):
        depth_all[i_time, i_lev, :] = -siglay[i_lev, :] * (h + zeta[i_time, :])

# 计算cell深度
print('计算cell深度')
depth_current = np.zeros((zeta.shape[0], siglay.shape[0], len(latc)), dtype=float)  # (time,lev,cell_num)
for i_lev in range(siglay.shape[0]):
    print(i_lev)
    for i_cell_subfix in range(len(lon)):
        depth_current[0, i_lev, i_cell_subfix] = np.average([
            depth_all[0, i_lev, nv[0, i_cell_subfix]],
            depth_all[0, i_lev, nv[1, i_cell_subfix]],
            depth_all[0, i_lev, nv[2, i_cell_subfix]]
        ])
        depth_current[:, i_lev, i_cell_subfix] = depth_current[0, i_lev, i_cell_subfix]

# 创建规则网格
print('创建网格')
lon_min, lon_max = 128.5, 132.5
lat_min, lat_max = 28, 32
grid_x, grid_y = np.mgrid[lon_min:lon_max:200j, lat_min:lat_max:200j]  # 提高分辨率到200j以获得更好的平滑效果

print('读取海岸线数据')
shp_path = r'D:\Documents\Data\海岸线数据\land_polygons.shp'
reader = shpreader.Reader(shp_path)

# 创建空间索引
geoms = list(reader.geometries())
tree = STRtree(geoms)

# 创建目标区域边界框
target_bbox = box(lon_min, lat_min, lon_max, lat_max)
prep_bbox = prep(target_bbox)  # 预处理加速查询

# 查询相交的几何体
intersect_indices = tree.query(target_bbox, predicate="intersects")
valid_geoms = [geoms[i] for i in intersect_indices]

# 使用更高效的方法计算陆地掩膜
print('使用matplotlib.Path计算陆地掩膜')
land_mask = np.zeros(grid_x.shape, dtype=bool)  # 初始化陆地掩膜为False

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

# 筛选数据（只取表面层）
print('筛选表面层数据')
lonc_new = []
latc_new = []
u_new = []
v_new = []
for i in range(len(latc)):
    if lon_min - 0.5 < lonc[i] < lon_max + 0.5 and lat_min - 0.5 < latc[i] < lat_max + 0.5:
        lonc_new.append(lonc[i])
        latc_new.append(latc[i])
        # 只取表面层（第一层）
        u_new.append(u[0, 0, i])  # [时间步, 层, 节点]
        v_new.append(v[0, 0, i])
lonc_new = np.array(lonc_new)
latc_new = np.array(latc_new)
u_new = np.array(u_new)
v_new = np.array(v_new)

# 处理时间数据
print('处理时间数据')
times_all = []
for i_time in range(len(times)):
    temp = ''.join(times[i_time, j].decode('utf-8') for j in range(len(times[0])))
    times_all.append((pd.to_datetime(temp) + timedelta(hours=8)).strftime('%Y-%m-%d %H时'))

# 设置投影和范围
print('设置投影与绘图范围')
proj = ccrs.PlateCarree()
extents = [lon_min, lon_max, lat_min, lat_max]

for i_time in range(u.shape[0]):
    print(f"正在处理时间步: {i_time + 1}/{u.shape[0]} - {times_all[i_time]}")

    # 创建图形 - 设置背景为灰色
    fig = plt.figure(figsize=(12, 10), dpi=150)

    ax = plt.axes(projection=proj)
    ax.set_extent(extents, crs=proj)
    ax.set_facecolor('lightgray')  # 设置坐标轴区域背景为灰色
    ax.set_title(times_all[i_time] + " 表面流速", fontsize=16)

    # 获取当前时间步的表面数据
    current_u = u_new  # 已经是表面数据
    current_v = v_new

    # 计算流速大小
    speed = np.sqrt(current_u ** 2 + current_v ** 2)

    # 插值到规则网格
    points_data = np.column_stack((lonc_new, latc_new))
    grid_speed = griddata(points_data, speed, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_u = griddata(points_data, current_u, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_v = griddata(points_data, current_v, (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 应用陆地掩膜
    grid_speed[land_mask] = np.nan
    grid_u[land_mask] = np.nan
    grid_v[land_mask] = np.nan

    # 计算单位向量（用于等长箭头）
    magnitude = np.sqrt(grid_u ** 2 + grid_v ** 2)
    grid_u_unit = np.where(magnitude > 0, grid_u / magnitude, 0)
    grid_v_unit = np.where(magnitude > 0, grid_v / magnitude, 0)

    # 使用自定义颜色映射 - 设置陆地掩膜区域为灰色
    ocean_cmap = plt.get_cmap('jet').copy()
    # ocean_cmap.set_bad('lightgray')  # 设置陆地区域使用灰色

    cf = ax.contourf(grid_x, grid_y, grid_speed, levels=50, cmap=ocean_cmap, transform=proj, cmin=0.,
                     cmax=np.nanmax(grid_speed))

    # 添加colorbar到右侧
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
    cbar.set_label('流速 (m/s)', fontsize=12)

    # 只绘制相关几何体 - 移除填充，保留边界线
    for geom in valid_geoms:
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            facecolor='none',  # 移除填充
            edgecolor='black',  # 保留黑色边界
            zorder=5
        )

    # 绘制等长流速箭头（只绘制在海洋区域）
    skip = 5  # 箭头密度
    # 创建海洋掩膜（非陆地）
    ocean_mask = ~land_mask
    # 创建索引数组
    i_indices = np.arange(grid_x.shape[0])[::skip]
    j_indices = np.arange(grid_x.shape[1])[::skip]

    # 收集有效点（海洋区域）
    valid_x = []
    valid_y = []
    valid_u = []
    valid_v = []

    for i in i_indices:
        for j in j_indices:
            if ocean_mask[i, j] and not np.isnan(grid_u_unit[i, j]) and not np.isnan(grid_v_unit[i, j]):
                valid_x.append(grid_x[i, j])
                valid_y.append(grid_y[i, j])
                valid_u.append(grid_u_unit[i, j])
                valid_v.append(grid_v_unit[i, j])

    if valid_x:  # 确保有有效点
        quiver = ax.quiver(
            valid_x,
            valid_y,
            valid_u,
            valid_v,
            scale=7.5,  # 控制箭头长度
            scale_units='inches',
            color='black',
            width=0.002,
            transform=proj,
            zorder=4  # 在陆地和海岸线之上
        )

    # 添加网格线
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # 保存图像
    filename = f"output/FVCOM/surface_flow_{i_time:02d}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.close()

    print(f"已保存: {filename}")

print("所有时间步处理完成")