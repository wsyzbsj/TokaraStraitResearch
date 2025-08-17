# 导入库
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

# 设置
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
lon_max = 130.2
lon_min = 129.6
lat_max = 30.2
lat_min = 29.85

# 创建自定义colormap - 海底使用棕色
brown_color = [0.55, 0.27, 0.07]  # RGB值表示深棕色
my_cmap = mcolors.ListedColormap([brown_color])

# 读取水深数据
gebco_file = r'D:\Documents\Data\gebco_2024\GEBCO_2024.nc'
depth_data = netCDF4.Dataset(gebco_file, 'r')
depth_lon = depth_data.variables['lon'][:]
depth_lat = depth_data.variables['lat'][:]
depth_elevation = depth_data.variables['elevation'][:]

# 获取目标区域数据 - 修正索引错误
lon_indices = np.where((depth_lon >= lon_min) & (depth_lon <= lon_max))[0]
lat_indices = np.where((depth_lat >= lat_min) & (depth_lat <= lat_max))[0]

# 确保有足够的点
if len(lon_indices) < 2 or len(lat_indices) < 2:
    raise ValueError("所选区域内数据点不足，请检查经纬度范围")

lon_used = depth_lon[lon_indices]
lat_used = depth_lat[lat_indices]
depth_used = depth_elevation[lat_indices[0]:lat_indices[-1]+1, lon_indices[0]:lon_indices[-1]+1]

# 下采样数据以提高性能
downsample_factor = 4  # 平衡细节和性能
X, Y = np.meshgrid(lon_used[::downsample_factor], lat_used[::downsample_factor])
Z = depth_used[::downsample_factor, ::downsample_factor]

# 创建3D图形
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# 设置垂直缩放因子（增强地形起伏）
vertical_exaggeration = 50

# 创建颜色数组 - 仅海底使用棕色
# 陆地（Z>=0）设置为完全透明（不显示）
# 海洋（Z<0）使用棕色
colors = np.zeros((Z.shape[0], Z.shape[1], 4))  # RGBA数组
land_mask = (Z >= 0)
ocean_mask = (Z < 0)

# 陆地设置为完全透明（不显示）
colors[land_mask] = [0.0, 0.0, 0.0, 0.0]  # RGBA: 全透明

# 海洋区域使用棕色
colors[ocean_mask] = [*brown_color, 1.0]  # RGBA: 棕色，完全不透明

# 绘制3D地形表面 - 修正底部移位
surf = ax.plot_surface(
    X, Y, Z * vertical_exaggeration,
    facecolors=colors,
    rstride=1, cstride=1,
    linewidth=0.1,  # 添加细线增强地形特征
    edgecolor='black',  # 黑色线条增强轮廓
    alpha=0.9,  # 轻微透明度增强深度感
    antialiased=True,  # 启用抗锯齿
    shade=True  # 启用阴影
)

# 设置视角 - 优化3D效果
ax.view_init(elev=38, azim=-60)  # 优化视角减少变形

# 设置坐标轴标签
ax.set_xlabel('经度 (°E)', labelpad=15, fontsize=10)
ax.set_ylabel('纬度 (°N)', labelpad=15, fontsize=10)
ax.set_zlabel('高程 (m)', labelpad=10, fontsize=10)

# 设置图形标题
ax.set_title('三维海底地形图', pad=15, fontsize=14, weight='bold')

# 精确设置坐标轴范围防止图形溢出
ax.set_xlim(np.min(X), np.max(X))
ax.set_ylim(np.min(Y), np.max(Y))

# 计算合理的Z轴范围
if np.any(ocean_mask):
    min_depth = np.min(Z[ocean_mask]) * vertical_exaggeration
else:
    min_depth = np.min(Z) * vertical_exaggeration
max_elevation = np.max(Z) * vertical_exaggeration if np.any(Z > 0) else 0

# 设置Z轴范围，确保底部不超出
ax.set_zlim(min_depth , 0)

# 调整布局增加边距
plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95)

# 添加网格线增强可读性
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.grid(True, linestyle='--', alpha=0.3)

# 保存高质量图形
plt.savefig('output/三维海底地形图.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
plt.show()