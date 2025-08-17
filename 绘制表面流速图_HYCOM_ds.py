import toml
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime, timedelta
import cartopy.io.shapereader as shpreader
from shapely.geometry import box, Point
from shapely.prepared import prep
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

# 读取流场数据
dataset = netCDF4.Dataset(r'D:\Documents\Data\HYCOM\uv3z_2024071709.nc4', 'r')
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
u = dataset.variables['water_u'][:]     # (time,depth,lat,lon)
v = dataset.variables['water_v'][:]     # (time,depth,lat,lon)
time = dataset.variables['time'][:][0]

# 修复1：确保坐标切片一致
print("筛选数据...")
lat_start = np.argwhere(lat > lat_min-0.5)[0][0]
lat_end = np.argwhere(lat < lat_max+0.5)[-1][0] + 1  # +1因为切片是左闭右开
lon_start = np.argwhere(lon > lon_min-0.5)[0][0]
lon_end = np.argwhere(lon < lon_max+0.5)[-1][0] + 1

lat_used = lat[lat_start:lat_end]
lon_used = lon[lon_start:lon_end]
u_used = u[0, 0, lat_start:lat_end, lon_start:lon_end]
v_used = v[0, 0, lat_start:lat_end, lon_start:lon_end]
dataset.close()

# # 修复NaN值
# u_used = np.nan_to_num(u_used, nan=0.0)
# v_used = np.nan_to_num(v_used, nan=0.0)

# 计算流速大小
print("计算流速...")
speed = np.sqrt(u_used ** 2 + v_used ** 2)

# 创建规则网格用于插值
print('创建网格...')
grid_x, grid_y = np.mgrid[lon_min:lon_max:200j, lat_min:lat_max:200j]

print('读取海岸线数据...')
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
print('计算陆地掩膜...')
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
print('设置投影与绘图范围...')
proj = ccrs.PlateCarree()
extents = [lon_min, lon_max, lat_min, lat_max]

time_str = (datetime(2000,1,1)+timedelta(hours=time)).strftime('%Y年%m月%d日%H时')

# 创建图形
fig = plt.figure(figsize=(12, 10), dpi=150)
ax = plt.axes(projection=proj)
ax.set_extent(extents, crs=proj)
ax.set_facecolor('lightgray')  # 设置坐标轴区域背景为灰色
ax.set_title(time_str + " HYCOM分析数据表面流速", fontsize=16)

# 修复2：创建原始数据网格点
lon_data, lat_data = np.meshgrid(lon_used, lat_used)
points_data = np.column_stack((lon_data.ravel(), lat_data.ravel()))

# 展平流速数据
speed_flat = speed.ravel()
u_flat = u_used.ravel()
v_flat = v_used.ravel()

# 插值到规则网格
print('插值数据到规则网格...')
grid_speed = griddata(points_data, speed_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)
grid_u = griddata(points_data, u_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)
grid_v = griddata(points_data, v_flat, (grid_x, grid_y), method='linear', fill_value=np.nan)

# 应用陆地掩膜
grid_speed[land_mask] = np.nan
grid_u[land_mask] = np.nan
grid_v[land_mask] = np.nan

# 计算单位向量（用于等长箭头）
magnitude = np.sqrt(grid_u ** 2 + grid_v ** 2)
grid_u_unit = np.where(magnitude > 0, grid_u / magnitude, 0)
grid_v_unit = np.where(magnitude > 0, grid_v / magnitude, 0)

# 使用自定义颜色映射
ocean_cmap = plt.get_cmap('jet').copy()
ocean_cmap.set_bad('lightgray')  # 设置陆地区域使用灰色

cf = ax.contourf(grid_x, grid_y, grid_speed, levels=50, cmap=ocean_cmap, transform=proj, cmin=0.,
                 cmax=np.nanmax(grid_speed))

# 添加colorbar到右侧
cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
cbar.set_label('流速 (m/s)', fontsize=12)

# 绘制海岸线
for geom in valid_geoms:
    ax.add_geometries(
        [geom],
        crs=ccrs.PlateCarree(),
        facecolor='none',  # 移除填充
        edgecolor='black',  # 保留黑色边界
        zorder=5
    )

# 绘制等长流速箭头（只绘制在海洋区域）
print('绘制流速箭头...')
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
filename = f"output/HYCOM/surface_flow_HYCOM_{time_str.replace(' ', '_').replace(':', '')}_表面流速.png"
plt.savefig(filename, bbox_inches='tight', dpi=150)
plt.close()

print(f"已保存: {filename}")