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
ax_main = fig.add_axes([0.05, 0.15, 0.55, 0.7], projection=ccrs.PlateCarree())
ax_main.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
ax_main.set_title('Full Domain with Zoom Area', fontsize=14)

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