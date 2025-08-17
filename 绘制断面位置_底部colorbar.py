import glob
import netCDF4
import numpy as np
from geopy.distance import geodesic as gd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator
import cartopy.crs as ccrs
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'

# 计算大圆距离
def get_distance_from_location(p1: tuple, p2: tuple):
    return gd(p1, p2).km

# 创建colorbar
# 定义颜色
color_start = np.array([0.8941, 0.7529, 0.4235])
color_mid = np.array([1, 1, 1])
color_end = np.array([0, 0.1294, 0.7098])
# 定义两段的数量
n1 = 8
n2 = 28
total_n = n1 + n2
pos_mid = n1 / total_n
# 构建segmentdata
cdict = {
    'red': [
        (0, color_start[0], color_start[0]),
        (pos_mid, color_mid[0], color_mid[0]),
        (1, color_end[0], color_end[0])
    ],
    'green': [
        (0, color_start[1], color_start[1]),
        (pos_mid, color_mid[1], color_mid[1]),
        (1, color_end[1], color_end[1])
    ],
    'blue': [
        (0, color_start[2], color_start[2]),
        (pos_mid, color_mid[2], color_mid[2]),
        (1, color_end[2], color_end[2])
    ]
}
# 创建colormap，总离散级别数设为total_n（也可以设大一些，比如256，但位置比例不变）
my_cmap = mcolors.LinearSegmentedColormap('my_cmap', cdict, N=total_n)

# 最大深度
max_y = 1200

# 设置断面起点和终点
point_start_loc = (28.5, 129.5)
point_end_loc = (31, 130.5)
distance = get_distance_from_location(point_start_loc, point_end_loc)
print(f"断面总长度: {distance:.2f} km")

# 插值出距离相等的点
num_points = 100
lat_list = np.linspace(point_start_loc[0], point_end_loc[0], num=num_points)
lon_list = np.linspace(point_start_loc[1], point_end_loc[1], num=num_points)
points = np.column_stack((lat_list, lon_list))

# 绘图
gebco_file = r'D:\Documents\Data\gebco_2024\GEBCO_2024.nc'
depth_data = netCDF4.Dataset(gebco_file, 'r')
depth_lon = depth_data.variables['lon'][:]
depth_lat = depth_data.variables['lat'][:]
depth_elevation = depth_data.variables['elevation'][:]

fig = plt.figure(figsize=(6, 10))
ax = fig.add_axes([0.1, 0.15, 0.8, 0.7], projection=ccrs.PlateCarree())  # 调整地图位置，为底部colorbar留空间

# 设置地图范围
min_lon = min(lon_list) - 1
max_lon = max(lon_list) + 1
min_lat = min(lat_list) - 1
max_lat = max(lat_list) + 1
ax.set_extent([min_lon, max_lon,min_lat, max_lat],crs=ccrs.PlateCarree())
lon_used = depth_lon[np.argwhere(depth_lon>min_lon)[0][0]:np.argwhere(depth_lon<max_lon)[-1][0]]
lat_used = depth_lat[np.argwhere(depth_lat>min_lat)[0][0]:np.argwhere(depth_lat<max_lat)[-1][0]]
depth_used = depth_elevation[np.argwhere(depth_lat>min_lat)[0][0]:np.argwhere(depth_lat<max_lat)[-1][0],
                            np.argwhere(depth_lon>min_lon)[0][0]:np.argwhere(depth_lon<max_lon)[-1][0]]


# 添加地图要素
im = ax.pcolormesh(lon_used, lat_used, -depth_used,
                  cmap=my_cmap,
                  vmin=0, vmax=2000,
                  shading='auto',  # 自动处理网格点
                  transform=ccrs.PlateCarree())  # 确保正确的坐标变换

# 创建与地图宽度对齐的底部 colorbar
# 计算地图轴的位置和尺寸
ax_pos = ax.get_position()
cax_width = ax_pos.width  # 与地图宽度相同
cax_height = 0.02  # colorbar 高度
cax_left = ax_pos.x0  # 与地图左侧对齐
cax_bottom = ax_pos.y0 - 0.05  # 在地图下方，留出5%的间距

# 创建 colorbar 轴
cax = fig.add_axes([cax_left, cax_bottom, cax_width, cax_height])
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')  # 水平方向

# 确保颜色范围设置正确
im.set_clim(0, 2000)

# 设置刻度
cbar.ax.xaxis.set_major_locator(MultipleLocator(500))  # 水平 colorbar 使用 xaxis
cbar.ax.xaxis.set_minor_locator(MultipleLocator(100))  # 水平 colorbar 使用 xaxis
cbar.ax.tick_params(axis='x', which='minor', labelbottom=False)  # 隐藏副刻度标签

# 添加单位（m）到右侧
cbar.ax.set_title("m", loc="right", fontsize=10, pad=5)  # 调整位置和间距

# 创建陆地掩码
land_mask = (depth_used >= 0.)  # 深度≥0表示陆地

# 创建陆地数据（1表示陆地，NaN表示海洋）
land_data = np.where(land_mask, 1, np.nan)

# 使用绿色填充陆地
land_cmap = mcolors.ListedColormap(['green'])
ax.pcolormesh(lon_used, lat_used, land_data,
              cmap=land_cmap,
              shading='auto',
              transform=ccrs.PlateCarree(),
              alpha=0.7)  # 设置透明度使地图更美观

# 添加经纬度
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linestyle='-.')

# 绘制断面线
ax.plot(lon_list, lat_list, 'r-', linewidth=2, transform=ccrs.Geodetic())
ax.contour(lon_used, lat_used, depth_used, linewidths=0.5, colors='black', levels=[0])

# 标记起点和终点
ax.plot(point_start_loc[1], point_start_loc[0], 'go', markersize=10,
        transform=ccrs.PlateCarree(), label='起点')
ax.plot(point_end_loc[1], point_end_loc[0], 'ro', markersize=10,
        transform=ccrs.PlateCarree(), label='终点')

# 添加标注
ax.text(point_start_loc[1] + 0.1, point_start_loc[0], f'起点\n{point_start_loc[0]}°N, {point_start_loc[1]}°E',
        transform=ccrs.PlateCarree(), fontsize=9)
ax.text(point_end_loc[1] + 0.1, point_end_loc[0], f'终点\n{point_end_loc[0]}°N, {point_end_loc[1]}°E',
        transform=ccrs.PlateCarree(), fontsize=9)

ax.set_title('流速截面位置', fontsize=14, pad=20) # pad设置间距
plt.savefig('output/流速截面图.png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.show()