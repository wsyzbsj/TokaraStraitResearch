import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Rectangle, ConnectionPatch

# 设置
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
lon_max = 154.
lon_min = 98.
lat_max = 53.
lat_min = -2.

# 更新后的局部放大区域设置
zoom_lon_min, zoom_lon_max = 128.5, 132.5
zoom_lat_min, zoom_lat_max = 28, 32

# 读取FVCOM输出文件
dataset = netCDF4.Dataset(r'D:\Documents\code\Pycharm\FVCOM_analysis\kuroshio_2024_0717\kuroshio_2024_0717.nc')
nv = dataset.variables['nv'][:]  # (3,441027) 网格点序号
lonc = dataset.variables['lonc'][:]  # (441027,)  网格中心
latc = dataset.variables['latc'][:]  # (441027,)  网格中心
lon = dataset.variables['lon'][:]  # (226031,)  网格点
lat = dataset.variables['lat'][:]  # (226031,)  网格点

# 创建图形和主图
fig = plt.figure(figsize=(16, 10), dpi=600)  # 调整为更宽的图形

# 主图设置 (左侧)
ax_main = fig.add_axes([0.05, 0.05, 0.40, 0.9], projection=ccrs.PlateCarree()) # 左底宽高
ax_main.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
# ax_main.set_title('Full Domain with Zoom Area', fontsize=14)

# 绘制所有网格
for i in range(lonc.shape[0]):
    point1_lat = lat[nv[0, i] - 1]
    point1_lon = lon[nv[0, i] - 1]
    point2_lat = lat[nv[1, i] - 1]
    point2_lon = lon[nv[1, i] - 1]
    point3_lat = lat[nv[2, i] - 1]
    point3_lon = lon[nv[2, i] - 1]
    lats = np.array([point1_lat, point2_lat, point3_lat, point1_lat])
    lons = np.array([point1_lon, point2_lon, point3_lon, point1_lon])
    ax_main.plot(lons, lats, color='black', linewidth=0.1)

# 在主图上添加框选区域
rect = Rectangle((zoom_lon_min, zoom_lat_min),
                 zoom_lon_max - zoom_lon_min,
                 zoom_lat_max - zoom_lat_min,
                 linewidth=1.5, edgecolor='red',
                 facecolor='none', zorder=10,
                 transform=ccrs.PlateCarree())
ax_main.add_patch(rect)

# 添加放大子图 (右侧)
ax_zoom = fig.add_axes([0.55, 0.05, 0.40, 0.9], projection=ccrs.PlateCarree())
ax_zoom.set_extent([zoom_lon_min, zoom_lon_max, zoom_lat_min, zoom_lat_max],crs=ccrs.PlateCarree())
# ax_zoom.set_title('Zoomed Area ('+str(zoom_lat_min)+'-'+str(zoom_lat_max)+'°N, '+str(zoom_lon_min)+'-'+str(zoom_lon_max)+'°E)', fontsize=14)

# 只绘制放大区域内的网格
for i in range(lonc.shape[0]):
    if (zoom_lon_min <= lonc[i] <= zoom_lon_max and
            zoom_lat_min <= latc[i] <= zoom_lat_max):
        point1_lat = lat[nv[0, i] - 1]
        point1_lon = lon[nv[0, i] - 1]
        point2_lat = lat[nv[1, i] - 1]
        point2_lon = lon[nv[1, i] - 1]
        point3_lat = lat[nv[2, i] - 1]
        point3_lon = lon[nv[2, i] - 1]

        lats = np.array([point1_lat, point2_lat, point3_lat, point1_lat])
        lons = np.array([point1_lon, point2_lon, point3_lon, point1_lon])

        ax_zoom.plot(lons, lats, color='black', linewidth=0.3)

# 添加地理特征
for ax in [ax_main, ax_zoom]:
    ax.coastlines('50m', linewidth=0.5)
    gl = ax.gridlines(draw_labels=True,linewidth=0.2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

# 添加连接线 - 正确连接左侧框角到右侧子图边框角
# 定义左侧框的四个角点（数据坐标）
main_corners = [
    (zoom_lon_min, zoom_lat_min),   # 左下角
    (zoom_lon_min, zoom_lat_max)    # 左上角
]

# 定义右侧子图的四个边框角（轴坐标）
zoom_corners = [
    (0, 0),  # 左下角
    (0, 1)  # 左上角
]

# 创建连接线
for i in range(2):
    con = ConnectionPatch(
        xyA=main_corners[i],  # 左侧框的角点（数据坐标）
        xyB=zoom_corners[i],  # 右侧子图的边框角（轴坐标）
        coordsA="data",  # A点使用数据坐标系
        coordsB="axes fraction",  # B点使用轴坐标系（0-1范围）
        axesA=ax_main,
        axesB=ax_zoom,
        color="red",
        linestyle="dashed",
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
plt.savefig('output/FVCOM网格.png', bbox_inches='tight', dpi=600)
plt.show()