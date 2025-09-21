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
temp = dataset.variables['temp'][:,0,:]  # 表面温度

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
    # ax.set_title(times_all[i_time] + " 表面温度与流速", fontsize=16)

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

    # 自定义colorbar
    ocean_cmap = plt.get_cmap('jet').copy()
    ocean_cmap.set_bad('lightgray')  # 设置陆地区域使用灰色

    # 绘制温度填色图 - 设置在最底层
    temp_cmap = plt.get_cmap('coolwarm')  # 使用适合温度的颜色映射
    cf_temp = ax.contourf(
        grid_x, grid_y, grid_temp,
        levels=np.linspace(27, 32, 50),
        cmap=ocean_cmap,
        transform=proj,
        zorder=1,  # 设置在最底层
        extend='both'
    )


    # 添加温度colorbar
    cbar_temp = plt.colorbar(cf_temp, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
    cbar_temp.set_label('温度 (℃)', fontsize=12)

    # 设置 colorbar 刻度
    tick_positions = np.arange(27, 32.1, 0.5)  # 从27到29，步长为0.5
    cbar_temp.set_ticks(tick_positions)

    # 只绘制相关几何体 - 移除填充，保留边界线
    for geom in valid_geoms:
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            facecolor='none',  # 移除填充
            edgecolor='black',  # 保留黑色边界
            zorder=5
        )

    # 绘制黑色流速箭头（调整箭头大小）
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
    valid_speed = []  # 用于过滤低速点

    for i in i_indices:
        for j in j_indices:
            if ocean_mask[i, j] and not np.isnan(grid_u[i, j]) and not np.isnan(grid_v[i, j]):
                # 过滤掉流速过小的点（小于0.05 m/s）
                if grid_speed[i, j] > 0.05:
                    valid_x.append(grid_x[i, j])
                    valid_y.append(grid_y[i, j])
                    valid_u.append(grid_u[i, j])
                    valid_v.append(grid_v[i, j])
                    valid_speed.append(grid_speed[i, j])

    quiver = None  # 初始化quiver变量
    if valid_x:  # 确保有有效点
        # 创建流速箭头图 - 使用黑色固定颜色
        # 调整箭头参数使箭头更明显
        quiver = ax.quiver(
            valid_x,
            valid_y,
            valid_u,
            valid_v,
            color='black',      # 统一黑色
            scale=5,            # 减小scale值使箭头变长
            scale_units='inches',
            width=0.003,        # 杆的宽度（稍微加粗）
            headwidth=3.0,       # 头部宽度
            headlength=4.0,      # 头部长度
            headaxislength=3.5,  # 头部轴长度
            transform=proj,
            zorder=4  # 在温度填色图之上，海岸线之下
        )

    # 添加网格线 - 设置在最底层
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', zorder=2)
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
    filename = f"output/FVCOM/surface_temp_flow_{i_time:02d}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.close()

    print(f"已保存: {filename}")

print("所有时间步处理完成")