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
from scipy.ndimage import gaussian_filter
import os

# 确保输出目录存在
os.makedirs("output/FVCOM", exist_ok=True)

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
dataset = netCDF4.Dataset(config['file_config']['fvcom_file'], 'r')
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
temp = dataset.variables['temp'][:,0,:]  # 表面温度

# 创建规则网格 - 保持较高分辨率但减少流线密度
print('创建网格')
lon_min, lon_max = 128.5, 132.5
lat_min, lat_max = 28, 32
grid_x, grid_y = np.mgrid[lon_min:lon_max:1000j, lat_min:lat_max:1000j]  # 250j分辨率

print('读取海岸线数据')
shp_path = config['file_config']['polygon_file']
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

# 分别获取目标区域内的节点索引和单元索引
print('获取目标区域内的节点和单元索引')
valid_node_indices = []  # 用于温度数据（节点数据）
valid_cell_indices = []  # 用于流速数据（单元数据）

# 获取有效节点索引（用于温度）
for i in range(len(lat)):
    if lon_min - 0.5 < lon[i] < lon_max + 0.5 and lat_min - 0.5 < lat[i] < lat_max + 0.5:
        valid_node_indices.append(i)
valid_node_indices = np.array(valid_node_indices)

# 获取有效单元索引（用于流速）
for i in range(len(latc)):
    if lon_min - 0.5 < lonc[i] < lon_max + 0.5 and lat_min - 0.5 < latc[i] < lat_max + 0.5:
        valid_cell_indices.append(i)
valid_cell_indices = np.array(valid_cell_indices)

# 处理时间数据
print('处理时间数据')
times_all = []
for i_time in range(len(times)):
    temp_str = ''.join(times[i_time, j].decode('utf-8') for j in range(len(times[0])))
    times_all.append((pd.to_datetime(temp_str) + timedelta(hours=8)).strftime('%Y-%m-%d %H时'))

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
    ax.set_title(times_all[i_time] + " 表面温度与流线", fontsize=16)

    # 获取当前时间步的数据
    current_u = u[i_time, 0, valid_cell_indices]  # 表面层u分量（单元数据）
    current_v = v[i_time, 0, valid_cell_indices]  # 表面层v分量（单元数据）
    current_temp = temp[i_time, valid_node_indices]  # 表面温度（节点数据）
    current_speed = np.sqrt(current_u ** 2 + current_v ** 2)  # 流速大小

    # 插值到规则网格
    # 温度插值（使用节点位置）
    points_temp = np.column_stack((lon[valid_node_indices], lat[valid_node_indices]))
    grid_temp = griddata(points_temp, current_temp, (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 流速插值（使用单元中心位置）
    points_flow = np.column_stack((lonc[valid_cell_indices], latc[valid_cell_indices]))
    grid_u = griddata(points_flow, current_u, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_v = griddata(points_flow, current_v, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_speed = griddata(points_flow, current_speed, (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 应用陆地掩膜
    grid_temp[land_mask] = np.nan
    grid_u[land_mask] = np.nan
    grid_v[land_mask] = np.nan
    grid_speed[land_mask] = np.nan

    # 平滑流速场 - 使用适中的平滑减少不连续性
    grid_u_smooth = gaussian_filter(grid_u, sigma=1.0)
    grid_v_smooth = gaussian_filter(grid_v, sigma=1.0)
    grid_speed_smooth = np.sqrt(grid_u_smooth**2 + grid_v_smooth**2)

    # 绘制温度填色图 - 设置在最底层
    temp_cmap = plt.get_cmap('coolwarm')  # 使用适合温度的颜色映射
    cf_temp = ax.contourf(
        grid_x, grid_y, grid_temp,
        levels=50,
        cmap=temp_cmap,
        transform=proj,
        zorder=1  # 设置在最底层
    )

    # 添加温度colorbar
    cbar_temp = plt.colorbar(cf_temp, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
    cbar_temp.set_label('温度 (℃)', fontsize=12)

    # 只绘制相关几何体 - 移除填充，保留边界线
    for geom in valid_geoms:
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            facecolor='lightgray',  # 移除填充
            edgecolor='black',  # 保留黑色边界
            zorder=10
        )

    # 创建海洋掩膜（非陆地）
    ocean_mask = ~land_mask

    # 创建海洋区域的流速网格（陆地区域设为np.nan）
    grid_u_ocean = np.where(ocean_mask, grid_u_smooth, np.nan)
    grid_v_ocean = np.where(ocean_mask, grid_v_smooth, np.nan)

    # 确保grid_speed_smooth中没有NaN值
    grid_speed_smooth_filled = np.nan_to_num(grid_speed_smooth, nan=np.nan)

    # 使用基于流速的加权采样生成起点 - 流速越快的地方起点越少，但不是完全没有
    num_start_points = 200  # 起点数量
    ocean_points = np.column_stack([grid_x[ocean_mask], grid_y[ocean_mask]])
    ocean_speeds = grid_speed_smooth_filled[ocean_mask]

    # 确保没有NaN值
    ocean_speeds = np.nan_to_num(ocean_speeds, nan=np.nan)

    # 计算权重：使用非线性函数，确保高速区域也有一定数量的起点
    max_speed = np.max(ocean_speeds) if np.max(ocean_speeds) > 0 else 1.0

    # 使用S形函数来平衡权重分布
    # 这个函数在低速区域权重较高，高速区域权重较低，但不会降到0
    def sigmoid_weight(speed, max_speed):
        # 将速度归一化到0-1范围
        normalized_speed = speed / max_speed
        # 使用S形函数的反函数，使高速区域权重降低但不会为0
        return 1.0 / (1.0 + np.exp(5.0 * (normalized_speed - 0.5)))

    weights = sigmoid_weight(ocean_speeds, max_speed)

    # 如果所有权重都是0，则使用均匀分布
    if np.sum(weights) == 0:
        weights = np.ones_like(weights)

    weights = weights / np.sum(weights)  # 归一化

    # 根据权重采样
    if len(ocean_points) > num_start_points:
        try:
            selected_indices = np.random.choice(
                len(ocean_points),
                size=num_start_points,
                replace=False,
                p=weights
            )
            start_points = ocean_points[selected_indices]
        except ValueError as e:
            print(f"权重采样错误: {e}, 使用均匀采样")
            selected_indices = np.random.choice(
                len(ocean_points),
                size=min(num_start_points, len(ocean_points)),
                replace=False
            )
            start_points = ocean_points[selected_indices]
    else:
        start_points = ocean_points

    # 绘制连续流线图
    if len(start_points) > 0:
        # 使用黑色流线，固定颜色
        stream = ax.streamplot(
            grid_x[:,0],  # x坐标（一维）
            grid_y[0,:],  # y坐标（一维）
            grid_u_ocean.T,  # 注意：需要转置以匹配streamplot的输入要求
            grid_v_ocean.T,
            color='black',  # 所有流线设为黑色
            linewidth=1.2,  # 线宽
            arrowsize=1.0,  # 箭头大小
            arrowstyle='->',  # 箭头样式
            minlength=0.1,  # 最小流线长度
            maxlength=100.0,  # 最大流线长度
            integration_direction='both',  # 双向积分
            start_points=start_points,  # 使用加权采样的起点
            broken_streamlines=False,    # 确保流线连续不断开
            density=1.2,  # 适中的密度值
            transform=proj,
            zorder=4
        )

    # 添加网格线 - 设置在最底层
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', zorder=2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # 保存图像
    filename = f"output/FVCOM/surface_temp_stream_{i_time:02d}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.close()

    print(f"已保存: {filename}")

print("所有时间步处理完成")