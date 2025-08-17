from datetime import datetime, timedelta
import numpy as np
import openpyxl
import ast
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geopy.distance import geodesic as gd
import os
import shutil
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

# 设置字体族，中文为微软雅黑，英文为Times New Roman
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
# 设置数学公式字体为stix
plt.rcParams['mathtext.fontset'] = 'stix'

# 读取HYCOM数据
dataset = netCDF4.Dataset(r'D:\Documents\Data\HYCOM\uv3z_2024.nc4', 'r')
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]

# 筛选区域内的网格
lon_max = 140
lon_min = 120
lat_max = 40
lat_min = 20
grid_step = 0.08

# 插值
# 创建目标网格
num_lon = int((lon_max - lon_min) / grid_step) + 1
num_lat = int((lat_max - lat_min) / grid_step) + 1
target_lon = np.linspace(lon_min, lon_max, num_lon)
target_lat = np.linspace(lat_min, lat_max, num_lat)
target_lon_grid, target_lat_grid = np.meshgrid(target_lon, target_lat)  # (num_lat, num_lon)

# 读取选定区域数据
lon_hycom = lon[np.argwhere(lon > lon_min)[0][0]:np.argwhere(lon < lon_max)[-1][0]]
lat_hycom = lat[np.argwhere(lat > lat_min)[0][0]:np.argwhere(lat < lat_max)[-1][0]]

u = dataset.variables['water_u'][:]
v = dataset.variables['water_v'][:]
depth = dataset.variables['depth'][:]

u_used = u[0, 0,
         np.argwhere(lat > lat_min)[0][0]:np.argwhere(lat < lat_max)[-1][0],
         np.argwhere(lon > lon_min)[0][0]:np.argwhere(lon < lon_max)[-1][0]]  # (深度,纬度,经度)

v_used = v[0, 0,
         np.argwhere(lat > lat_min)[0][0]:np.argwhere(lat < lat_max)[-1][0],
         np.argwhere(lon > lon_min)[0][0]:np.argwhere(lon < lon_max)[-1][0]]  # (深度,纬度,经度)

speed_used = np.sqrt(u_used ** 2 + v_used ** 2)  # 计算速度

# 读取水深数据
gebco_file = r'D:\Documents\Data\gebco_2024\GEBCO_2024.nc'
depth_data = netCDF4.Dataset(gebco_file, 'r')
depth_lon = depth_data.variables['lon'][:]
depth_lat = depth_data.variables['lat'][:]
depth_elevation = depth_data.variables['elevation'][:]

lon_used = depth_lon[np.argwhere(depth_lon > lon_min)[0][0]:np.argwhere(depth_lon < lon_max)[-1][0]]
lat_used = depth_lat[np.argwhere(depth_lat > lat_min)[0][0]:np.argwhere(depth_lat < lat_max)[-1][0]]
depth_used = depth_elevation[np.argwhere(depth_lat > lat_min)[0][0]:np.argwhere(depth_lat < lat_max)[-1][0],
             np.argwhere(depth_lon > lon_min)[0][0]:np.argwhere(depth_lon < lon_max)[-1][0]]

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

# 创建colormap
my_cmap = mcolors.LinearSegmentedColormap('my_cmap', cdict, N=total_n)

# 绘图
fig = plt.figure(figsize=(12, 6))
ax = fig.add_axes([0.1, 0.15, 0.8, 0.7], projection=ccrs.PlateCarree())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# 添加水深图
im = ax.pcolormesh(lon_used, lat_used, -depth_used,
                   cmap=my_cmap,
                   vmin=0, vmax=2000,
                   shading='auto',
                   transform=ccrs.PlateCarree())

# 创建与地图宽度对齐的底部 colorbar
ax_pos = ax.get_position()
cax_width = ax_pos.width
cax_height = 0.02
cax_left = ax_pos.x0
cax_bottom = ax_pos.y0 - 0.05

# 创建 colorbar 轴
cax = fig.add_axes([cax_left, cax_bottom, cax_width, cax_height])
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')

# 确保颜色范围设置正确
im.set_clim(0, 2000)

# 设置刻度
cbar.ax.xaxis.set_major_locator(MultipleLocator(500))
cbar.ax.xaxis.set_minor_locator(MultipleLocator(100))
cbar.ax.tick_params(axis='x', which='minor', labelbottom=False)

# 添加单位（m）到右侧
cbar.ax.set_title("m", loc="right", fontsize=10, pad=5)

# 创建陆地掩码
land_mask = (depth_used >= 0.0)  # 深度≥0表示陆地

# 创建陆地数据（1表示陆地，NaN表示海洋）
land_data = np.where(land_mask, 1, np.nan)

# 使用绿色填充陆地
land_cmap = mcolors.ListedColormap(['green'])
ax.pcolormesh(lon_used, lat_used, land_data,
              cmap=land_cmap,
              shading='auto',
              transform=ccrs.PlateCarree(),
              alpha=0.7)

# ====== 添加0.8 m/s等值线 ======
# 绘制0.8 m/s的等值线（红色实线）
contour = ax.contour(
    lon_hycom, lat_hycom, speed_used,
    levels=[0.8],  # 指定0.8 m/s的等值线
    colors='red',  # 红色线条
    linewidths=1.5,  # 线宽
    linestyles='solid',  # 实线
    transform=ccrs.PlateCarree()
)

# 添加等值线标签（显示"0.8 m/s"）
ax.clabel(
    contour,
    inline=True,  # 标签嵌入在线条中
    fontsize=10,  # 字体大小
    fmt='%1.1f m/s',  # 标签格式
    colors='red'  # 标签颜色
)
# =============================

ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
ax.set_title('黑潮流轴')
plt.savefig('output/黑潮流轴_Speed.png', dpi=600, bbox_inches='tight')
plt.close()