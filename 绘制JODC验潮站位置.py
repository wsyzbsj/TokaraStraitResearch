import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle, ConnectionPatch
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator,FixedLocator
import toml

# 创建colorbar
# 定义颜色
color_start = np.array([1, 1, 1])
color_end = np.array([12./255, 87./255, 160./255])#0c57a0
# 定义总数量
total_n = 20
# 构建segmentdata
cdict = {
    'red': [
        (0, color_start[0], color_start[0]),
        (1, color_end[0], color_end[0])
    ],
    'green': [
        (0, color_start[1], color_start[1]),
        (1, color_end[1], color_end[1])
    ],
    'blue': [
        (0, color_start[2], color_start[2]),
        (1, color_end[2], color_end[2])
    ]
}
# 创建colormap，总离散级别数设为total_n（也可以设大一些，比如256，但位置比例不变）
my_cmap = mcolors.LinearSegmentedColormap('my_cmap', cdict, N=total_n)

# 设置
config = toml.load('config/config.toml')
plt.rcParams['font.family'] = 'Times New Roman, SimHei'
plt.rcParams['mathtext.fontset'] = 'stix'
# 大区域经纬度范围
lat_min = config['zone_config']['lat_min']
lat_max = config['zone_config']['lat_max']
lon_min = config['zone_config']['lon_min']
lon_max = config['zone_config']['lon_max']
# 局部放大区域
zoom_lon_min = config['zone_config']['lon_min_zoom']
zoom_lon_max = config['zone_config']['lon_max_zoom']
zoom_lat_min = config['zone_config']['lat_min_zoom']
zoom_lat_max = config['zone_config']['lat_max_zoom']


# 读取水深数据文件
dataset = netCDF4.Dataset('/home/yzbsj/Data/其他数据/GEBCO/gebco_2024/GEBCO_2024.nc')
depth_lon = dataset.variables['lon'][:]
depth_lat = dataset.variables['lat'][:]
depth_elevation = -dataset.variables['elevation'][:]

# 筛选全部区域
lon_used = depth_lon[np.argwhere(depth_lon > lon_min)[0][0]:np.argwhere(depth_lon < lon_max)[-1][0]]
lat_used = depth_lat[np.argwhere(depth_lat > lat_min)[0][0]:np.argwhere(depth_lat < lat_max)[-1][0]]
depth_used = depth_elevation[np.argwhere(depth_lat > lat_min)[0][0]:np.argwhere(depth_lat < lat_max)[-1][0],
np.argwhere(depth_lon > lon_min)[0][0]:np.argwhere(depth_lon < lon_max)[-1][0]]

# 筛选放大区域
lon_zoom = depth_lon[np.argwhere(depth_lon > zoom_lon_min)[0][0]:np.argwhere(depth_lon < zoom_lon_max)[-1][0]]
lat_zoom = depth_lat[np.argwhere(depth_lat > zoom_lat_min)[0][0]:np.argwhere(depth_lat < zoom_lat_max)[-1][0]]
depth_zoom = depth_elevation[np.argwhere(depth_lat > zoom_lat_min)[0][0]:np.argwhere(depth_lat < zoom_lat_max)[-1][0],
np.argwhere(depth_lon > zoom_lon_min)[0][0]:np.argwhere(depth_lon < zoom_lon_max)[-1][0]]

# 创建图形和主图
fig = plt.figure(figsize=(9, 12), dpi=600)

# 主图设置 (上部)
ax_main = fig.add_axes([0.06, 0.70, 0.80, 0.25], projection=ccrs.PlateCarree()) # 左底宽高
ax_main.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
# ax_main.set_title('Full Domain with Zoom Area', fontsize=14)
# 绘制水深
im = ax_main.contourf(lon_used, lat_used, depth_used, levels=np.linspace(0,10000,21), cmap=my_cmap, vmin=0, vmax=10000)
# 添加框选区域
rect = Rectangle((zoom_lon_min, zoom_lat_min),
                 zoom_lon_max - zoom_lon_min,
                 zoom_lat_max - zoom_lat_min,
                 linewidth=1.5, edgecolor='black',
                 facecolor='none', zorder=10,
                 transform=ccrs.PlateCarree())
ax_main.add_patch(rect)

# 添加字体
ax_main.text(120, 30.5, 'China', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')               # 中国
ax_main.text(123, 28.5, 'East China Sea', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')    # 东中国海
ax_main.text(131, 33, 'Japan', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')               # 日本
ax_main.text(126, 33.5, 'Korea', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')               # 韩国
ax_main.text(135, 26, 'North Pacific', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')       # 东中国海
# 添加colorbar
# 创建 colorbar 轴
cax = fig.add_axes([0.87, 0.70, 0.05, 0.25])
cbar = fig.colorbar(im, cax=cax, orientation='vertical', label="深度 (m)")  # 垂直方向
# 确保颜色范围设置正确
im.set_clim(0, 10000)
# 设置刻度
cbar.ax.yaxis.set_major_locator(MultipleLocator(2000))  # 垂直 colorbar 使用 yaxis
cbar.ax.yaxis.set_minor_locator(MultipleLocator(500))  # 垂直 colorbar 使用 yaxis
cbar.ax.tick_params(axis='y', which='minor')  # 隐藏副刻度标签
cbar.ax.invert_yaxis()
# 创建陆地掩码
land_mask = (depth_used <= 0.)  # 深度≥0表示陆地
# 创建陆地数据（1表示陆地，NaN表示海洋）
land_data = np.where(land_mask, 1, np.nan)
# 使用浅灰色色填充陆地
land_cmap = mcolors.ListedColormap(['lightgray'])
ax_main.pcolormesh(lon_used, lat_used, land_data,
                   cmap=land_cmap,
                   shading='auto',
                   transform=ccrs.PlateCarree(),
                   alpha=1)  # 设置透明度使地图更美观
ax_main.contour(lon_used, lat_used, depth_used, levels=[0], colors='black', linewidths=0.5)

# 添加放大子图 (下部)
ax_zoom = fig.add_axes([0.06, 0.05, 0.80, 0.6], projection=ccrs.PlateCarree())
ax_zoom.set_extent([zoom_lon_min, zoom_lon_max, zoom_lat_min, zoom_lat_max], crs=ccrs.PlateCarree())
# 使用pcolormesh绘制水深
im_zoom = ax_zoom.pcolormesh(lon_zoom, lat_zoom, depth_zoom,
                             vmin=0, vmax=2000,
                             cmap=my_cmap,
                             shading='auto')
# 设置超过2000m的颜色为深蓝色
im_zoom.cmap.set_over('#1f3c73')
# 创建 colorbar 轴
cax = fig.add_axes([0.87, 0.05, 0.05, 0.6])
cbar = fig.colorbar(im_zoom, cax=cax, orientation='vertical', label='深度 (m)')
# 反转colorbar并设置刻度
cbar.ax.invert_yaxis()
# 设置主刻度和标签位置
major_ticks = [0, 500, 1000, 1500, 2000]
cbar.ax.yaxis.set_major_locator(FixedLocator(major_ticks))
# 设置副刻度位置
minor_ticks = np.arange(0, 2001, 100)
cbar.ax.yaxis.set_minor_locator(FixedLocator(minor_ticks))
# 隐藏副刻度标签
cbar.ax.tick_params(axis='y', which='minor')
# 创建陆地掩码和绘制
land_mask = (depth_zoom <= 0.)
land_data = np.where(land_mask, 1, np.nan)
land_cmap = mcolors.ListedColormap(['lightgray'])
ax_zoom.pcolormesh(lon_zoom, lat_zoom, land_data,
                   cmap=land_cmap,
                   shading='auto',
                   transform=ccrs.PlateCarree(),
                   alpha=1)
ax_zoom.contour(lon_zoom, lat_zoom, depth_zoom, levels=[0], colors='black', linewidths=0.5)

# 标记起点和终点
ax_zoom.plot(130.689, 31.024, 'o', markersize=10, transform=ccrs.PlateCarree(), label='HD20')
ax_zoom.plot(129.498, 28.378, 'o', markersize=10, transform=ccrs.PlateCarree(), label='HD22')
ax_zoom.plot(130.995, 30.732, 'o', markersize=10, transform=ccrs.PlateCarree(), label='HD21')
ax_zoom.plot(129.848, 29.842, 'o', markersize=10, transform=ccrs.PlateCarree(), label='HD28')
ax_zoom.plot(129.533, 28.317, 'o', markersize=10, transform=ccrs.PlateCarree(), label='MA73')
ax_zoom.plot(130.97, 30.467, 'o', markersize=10, transform=ccrs.PlateCarree(), label='MA75')
ax_zoom.text(130.641, 31.024, 'ODOMARI', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='right')
ax_zoom.text(129.45, 28.45, 'NAZE', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='right')
ax_zoom.text(130.946, 30.732, 'NISHINOOMOTE', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='right')
ax_zoom.text(129.8, 29.842, 'NAKANOSHIMA', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='right')
ax_zoom.text(129.6, 28.30, 'AMAMI', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='left')
ax_zoom.text(130.92, 30.467, 'TANEGASHIMA', transform=ccrs.PlateCarree(), fontsize=12, horizontalalignment='right')

# 添加标注
ax_zoom.text(130.7, 31.8, 'Kyushu', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')                           # 九州
ax_zoom.text(129.5, 28.2, 'Amamioshima', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')                      # 奄美大岛
ax_zoom.text(129.8, 29.5, 'Tokara Islands', transform=ccrs.PlateCarree(), fontsize=15, rotation=55, fontweight='bold')       # 吐噶喇群岛
ax_zoom.text(128.5, 30.3, 'East China Sea', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')                   # 东中国海
ax_zoom.text(130.7, 30.2, 'Yakushima', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')                   # 东中国海
ax_zoom.text(131.1, 30.6, 'Tanegashima', transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold')                   # 东中国海

# 添加地理特征
for ax in [ax_main, ax_zoom]:
    gl = ax.gridlines(draw_labels=True,linewidth=0.5, color='gray', alpha=1, linestyle='--')
    gl.bottom_labels = False
    gl.right_labels = False
    if ax == ax_main:
        gl.ylocator = FixedLocator([20, 23, 26, 29, 32, 35])

# 添加连接线 - 正确连接左侧框角到右侧子图边框角
# 定义左侧框的四个角点（数据坐标）
main_corners = [
    (zoom_lon_max, zoom_lat_max),   # 右上角
    (zoom_lon_min, zoom_lat_max)    # 左上角
]

# 定义右侧子图的四个边框角（轴坐标）
zoom_corners = [
    (1, 1),  # 右上角
    (0, 1)  # 左上角
]

# 创建连接线
for i in range(len(main_corners)):
    con = ConnectionPatch(
        xyA=main_corners[i],  # 左侧框的角点（数据坐标）
        xyB=zoom_corners[i],  # 右侧子图的边框角（轴坐标）
        coordsA="data",  # A点使用数据坐标系
        coordsB="axes fraction",  # B点使用轴坐标系（0-1范围）
        axesA=ax_main,
        axesB=ax_zoom,
        color="black",
        linestyle="-",
        linewidth=1.5,
        alpha=0.9,
        zorder=15  # 确保在最上层
    )
    fig.add_artist(con)

# 添加经纬度标签
# ax_main.set_xlabel('Longitude (°E)', fontsize=12)
# ax_main.set_ylabel('Latitude (°N)', fontsize=12)
# ax_zoom.set_xlabel('Longitude (°E)', fontsize=12)

# 调整布局并保存
plt.savefig('output/验潮站位置.png', dpi=600)
plt.show()