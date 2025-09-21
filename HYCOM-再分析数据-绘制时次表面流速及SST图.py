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
from matplotlib.patches import Rectangle

# 设置
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
config = toml.load('config/config.toml')
# 局部放大区域
lon_min = config['zone_config']['lon_min_zoom']
lon_max = config['zone_config']['lon_max_zoom']
lat_min = config['zone_config']['lat_min_zoom']
lat_max = config['zone_config']['lat_max_zoom']

# 读取流场数据
dataset = netCDF4.Dataset(r'D:\Documents\Data\HYCOM\uv3z_2024071709.nc4', 'r')
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
u = dataset.variables['water_u'][:]     # (time,depth,lat,lon)
v = dataset.variables['water_v'][:]     # (time,depth,lat,lon)
time = dataset.variables['time'][:][0]
# 筛选数据
lat_idx_min = np.argwhere(lat > lat_min-1)[0][0]
lat_idx_max = np.argwhere(lat < lat_max+1)[-1][0]
lon_idx_min = np.argwhere(lon > lon_min-1)[0][0]
lon_idx_max = np.argwhere(lon < lon_max+1)[-1][0]

lat_used = lat[lat_idx_min:lat_idx_max]
lon_used = lon[lon_idx_min:lon_idx_max]
u_used = u[0, 0, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]
v_used = v[0, 0, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]
dataset.close()

# 读取SST数据
dataset_sst = netCDF4.Dataset(r'D:\Documents\Data\HYCOM\ts3z_20240717.nc4', 'r')
lat_sst = dataset_sst.variables['lat'][:]
lon_sst = dataset_sst.variables['lon'][:]
sst = dataset_sst.variables['water_temp'][:]
# 筛选数据
lat_sst_idx_min = np.argwhere(lat_sst > lat_min-1)[0][0]
lat_sst_idx_max = np.argwhere(lat_sst < lat_max+1)[-1][0]
lon_sst_idx_min = np.argwhere(lon_sst > lon_min-1)[0][0]
lon_sst_idx_max = np.argwhere(lon_sst < lon_max+1)[-1][0]

lat_sst_used = lat_sst[lat_sst_idx_min:lat_sst_idx_max]
lon_sst_used = lon_sst[lon_sst_idx_min:lon_sst_idx_max]
sst_used = sst[0, 0, lat_sst_idx_min:lat_sst_idx_max, lon_sst_idx_min:lon_sst_idx_max]
dataset_sst.close()

# 创建规则网格
print('创建网格')
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

# 设置投影和范围
print('设置投影与绘图范围')
proj = ccrs.PlateCarree()
extents = [lon_min, lon_max, lat_min, lat_max]

time_str = (datetime(2000,1,1)+timedelta(hours=time)).strftime('%Y年%m月%d日%H时')

# 创建图形 - 设置背景为灰色
fig = plt.figure(figsize=(12, 10), dpi=150)

ax = plt.axes(projection=proj)
ax.set_extent(extents, crs=proj)
ax.set_facecolor('lightgray')  # 设置坐标轴区域背景为灰色
# ax.set_title(time_str + "HYCOM分析数据表面流速", fontsize=16)

# 获取当前时间步的表面数据
current_u = u_used  # 已经是表面数据
current_v = v_used

# 计算流速大小 (这里可能有NaN值，所以会有警告，但可以忽略)
speed = np.sqrt(current_u ** 2 + current_v ** 2)

# 创建流场数据的网格点
lon_grid, lat_grid = np.meshgrid(lon_used, lat_used)
points_data = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

# 创建SST数据的网格点
lon_sst_grid, lat_sst_grid = np.meshgrid(lon_sst_used, lat_sst_used)
points_sst_data = np.column_stack((lon_sst_grid.ravel(), lat_sst_grid.ravel()))

# 展平流场和SST数据
u_flat = current_u.ravel()
v_flat = current_v.ravel()
sst_flat = sst_used.ravel()

# 插值到规则网格
grid_sst = griddata(points_sst_data, sst_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)
grid_u = griddata(points_data, u_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)
grid_v = griddata(points_data, v_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)
for i in range(grid_sst.shape[0]):
    for j in range(grid_sst.shape[1]):
        if grid_sst[i, j] < 0:
            grid_sst[i, j] = np.nan
for i in range(grid_u.shape[0]):
    for j in range(grid_v.shape[1]):
        if (grid_u[i,j]**2+grid_v[i,j]**2)**0.5 > 2:
            grid_u[i,j] = np.nan
            grid_v[i,j] = np.nan

# 应用陆地掩膜
grid_u[land_mask] = np.nan
grid_v[land_mask] = np.nan
grid_sst[land_mask] = np.nan

# 使用自定义颜色映射
ocean_cmap = plt.get_cmap('jet').copy()
ocean_cmap.set_bad('lightgray')  # 设置陆地区域使用灰色

cf = ax.contourf(grid_x, grid_y, grid_sst, levels=np.linspace(27, 32, 50), cmap=ocean_cmap, transform=proj, extend='both')

# 添加colorbar到右侧
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
cbar.set_label('SST (℃)', fontsize=12)

# 设置 colorbar 刻度
tick_positions = np.arange(27, 32.1, 0.5)  # 从27到29，步长为0.5
cbar.set_ticks(tick_positions)

# 如果需要，可以格式化刻度标签
cbar.set_ticklabels([f'{tick:.2f}' for tick in tick_positions])

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
        if ocean_mask[i, j] and not np.isnan(grid_u[i, j]) and not np.isnan(grid_v[i, j]):
            valid_x.append(grid_x[i, j])
            valid_y.append(grid_y[i, j])
            valid_u.append(grid_u[i, j])
            valid_v.append(grid_v[i, j])

quiver = None
if valid_x:  # 确保有有效点
    quiver = ax.quiver(
        valid_x,
        valid_y,
        valid_u,
        valid_v,
        scale=5.0,  # 控制箭头长度
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

# 添加图例（quiverkey）在最上层
if quiver is not None:
    # 添加空白背景框
    bg_patch = Rectangle(
        (0.88, 0.0), 0.12, 0.05,  # (x, y), width, height
        transform=ax.transAxes,
        facecolor='white',
        alpha=0.8,  # 半透明
        edgecolor='black',
        linewidth=0.5,
        zorder=10  # 设置在最上层
    )
    ax.add_patch(bg_patch)

    # 添加quiverkey图例（显示1m/s向东流的箭头）
    ax.quiverkey(
        quiver,
        X=0.93, Y=0.025,    # 宽度的一半,高度的一半
        U=1,
        label='1 m/s',
        labelpos='E',
        coordinates='axes',
        fontproperties={'size': 10, 'weight': 'bold'},
        labelsep=0.1,  # 标签与箭头之间的距离
        color='black',
        zorder=11  # 设置在最上层
    )

# 保存图像
filename = f"output/HYCOM/surface_flow_HYCOM"+time_str+"表面流速.png"
plt.savefig(filename, bbox_inches='tight', dpi=150)
plt.close()

print(f"已保存: {filename}")