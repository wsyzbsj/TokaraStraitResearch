import toml
from geopy.distance import geodesic as gd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import cmocean
plt.rcParams['font.family'] = 'Times New Roman, SimHei'
plt.rcParams['mathtext.fontset'] = 'stix'

def haversine_vectorized(lat1, lon1, lat2, lon2): # 向量化的Haversine距离计算
    # 将十进制度数转化为弧度
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Haversine公式
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    r = 6371  # 地球平均半径，单位公里
    return c * r

if __name__ == "__main__":
    # 读取设置——注:用不到大区域故将lat_max等变量名给数据范围
    config = toml.load('config/config.toml')
    current_section = config['line_config']['current_section']['section']
    line_start_lat = config['line_config'][current_section]['start_lat']
    line_start_lon = config['line_config'][current_section]['start_lon']
    line_end_lat = config['line_config'][current_section]['end_lat']
    line_end_lon = config['line_config'][current_section]['end_lon']
    point_start_loc = (line_start_lat, line_start_lon)
    point_end_loc = (line_end_lat, line_end_lon)
    lat_max = max(point_start_loc[0], point_end_loc[0])+0.5
    lat_min = min(point_start_loc[0], point_end_loc[0])-0.5
    lon_max = max(point_start_loc[1], point_end_loc[1])+0.5
    lon_min = min(point_start_loc[1], point_end_loc[1])-0.5
    distance = gd(point_start_loc, point_end_loc).km
    print(f"断面总长度: {distance:.2f} km")

    # 插值出距离相等的点
    num_points = 100
    lat_list = np.linspace(point_start_loc[0], point_end_loc[0], num=num_points)
    lon_list = np.linspace(point_start_loc[1], point_end_loc[1], num=num_points)
    points = np.column_stack((lat_list, lon_list))

    if current_section=='TK' or current_section=='TokaraStrait':
        plot_x = lat_list
        x_label = 'Latitude (°N)'
    elif current_section=='PN':
        plot_x = lon_list
        x_label = 'Longitude (°E)'

    # === 1. 处理水深数据 ===
    depth_file = '/home/yzbsj/Data/海洋数据/HYCOM/expt_53.x_meanstd/regional_depth_11.nc'
    depth_data = netCDF4.Dataset(depth_file, 'r')
    depth_lon = depth_data.variables['Longitude'][:]
    depth_lat = depth_data.variables['Latitude'][:]
    depth_elevation = depth_data.variables['depth'][0,:]
    depth_data.close()
    # 处理经纬度信息
    loc = []
    for i in range(depth_lon.shape[0]):
        loc.append([])
        for j in range(depth_lon.shape[1]):
            loc[i].append((depth_lat[i,j],depth_lon[i, j])) # 纬度,经度

    # 计算线段上的深度
    seafloor_depth= []
    for i_point in range(len(points)):
        # 假设您有二维的经纬度数组（而不是列表的列表）
        # 如果您的数据是列表的列表，可以转换为NumPy数组
        depth_lat_array = np.array([[point[0] for point in row] for row in loc])
        depth_lon_array = np.array([[point[1] for point in row] for row in loc])

        # 目标点
        target_lat, target_lon = points[i_point,0],points[i_point,1]

        # 计算所有点与目标点的距离
        distances = haversine_vectorized(target_lat, target_lon, depth_lat_array, depth_lon_array)

        # 找到最小距离的索引
        min_index = np.unravel_index(np.argmin(distances), distances.shape)
        seafloor_depth.append(depth_elevation[min_index[0],min_index[1]])

    seafloor_depth = np.array(seafloor_depth)

    # 首先获取有效数据的索引和值
    valid_indices = np.where(~np.isnan(seafloor_depth))[0]  # 获取未屏蔽的索引
    valid_values = seafloor_depth[valid_indices]  # 获取未屏蔽的值

    # 获取所有需要插值的索引（被屏蔽的索引）
    masked_indices = np.where(np.isnan(seafloor_depth))[0]

    # 使用 np.interp 进行线性插值
    interpolated_values = np.interp(
        masked_indices,           # 需要插值的x坐标（索引位置）
        valid_indices,            # 已知数据的x坐标（有效数据的索引位置）
        valid_values              # 已知数据的y坐标（有效数据的值）
    )

    # 创建一个新的数组，将插值结果填充到被屏蔽的位置
    result = seafloor_depth.copy()  # 创建原始数组的副本
    result[masked_indices] = interpolated_values  # 用插值结果替换被屏蔽的值
    seafloor_depth = result

    # 设置绘图最深深度
    # max_y = (seafloor_depth.max()//100+1)*100 # 最大深度
    max_y = 1000 # 最大深度

    # === 2. 读取数据 ===
    spring_data = np.load('output/HYCOM/'+current_section+'/ave_spring_velocity_section.npz')
    summer_data = np.load('output/HYCOM/'+current_section+'/ave_summer_velocity_section.npz')
    autumn_data = np.load('output/HYCOM/'+current_section+'/ave_autumn_velocity_section.npz')
    winter_data = np.load('output/HYCOM/'+current_section+'/ave_winter_velocity_section.npz')

    # 设置最大最小值
    max_y_loc = np.min([np.argwhere(spring_data['y']<=max_y)[-1,0],len(spring_data['y'])-2])
    vmax = np.ceil(np.max([np.nanmax(spring_data['data'][:,:max_y_loc+2]), np.nanmax(summer_data['data'][:,:max_y_loc+2]),
                           np.nanmax(autumn_data['data'][:,:max_y_loc+2]), np.nanmax(winter_data['data'][:,:max_y_loc+2])])*10)/10
    vmin = np.floor(np.min([np.nanmin(spring_data['data'][:,:max_y_loc+2]), np.nanmin(summer_data['data'][:,:max_y_loc+2]),
                            np.nanmin(autumn_data['data'][:,:max_y_loc+2]), np.nanmin(winter_data['data'][:,:max_y_loc+2])])*10)/10

    # 创建图形并调整布局，为底部colorbar腾出空间
    fig = plt.figure(figsize=(10, 6))
    plt.subplots_adjust(bottom=0.08, top=0.90, right=0.85)  # 调整底部边距
    plt.suptitle('Water normal velocity at '+current_section+' section', fontsize=20)

    # 春季
    ax1 = plt.subplot(221)
    # 创建正确的二维网格坐标
    X_grid, Y_grid = np.meshgrid(plot_x, spring_data['y'][:max_y_loc+2], indexing='ij')  # 确保二维网格

    # 创建颜色映射
    cmap = cmocean.cm.balance
    cmap.set_bad('lightgray', 1.0)

    # 绘制水温数据 (使用二维坐标)
    mesh = ax1.contourf(
        X_grid, Y_grid, spring_data['data'][:,:max_y_loc+2],  # 使用二维坐标和流速数据
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        levels=np.arange(vmin, vmax+0.05, 0.05),
        zorder=0
    )

    # 添加等值线
    ax1.contour(
        X_grid, Y_grid, spring_data['data'][:,:max_y_loc+2],
        levels=[0],
        colors='k',
        linewidths=0.5,
        linestyles='--',
        zorder=0
    )

    # 绘制海床地形 (使用一维距离和深度)
    ax1.plot(plot_x, seafloor_depth, 'k-', linewidth=2, label='海床地形')
    ax1.fill_between(plot_x, seafloor_depth, max_y, color='lightgray', zorder=10)

    # 添加海平面线
    ax1.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

    # 设置坐标轴
    ax1.invert_yaxis()  # 深度向下增加
    ax1.set_ylim(max_y, 0)
    ax1.set_ylabel('Depth (m)', fontsize=12)
    ax1.set_title('Spring', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.75)

    # 夏季
    ax2 = plt.subplot(222)
    # 创建正确的二维网格坐标
    X_grid, Y_grid = np.meshgrid(plot_x, summer_data['y'][:max_y_loc+2], indexing='ij')  # 确保二维网格

    # 绘制水温数据 (使用二维坐标)
    mesh = ax2.contourf(
        X_grid, Y_grid, summer_data['data'][:,:max_y_loc+2],  # 使用二维坐标和流速数据
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        levels=np.arange(vmin, vmax+0.05, 0.05),
        zorder=0
    )

    # 添加等值线
    ax2.contour(
        X_grid, Y_grid, summer_data['data'][:,:max_y_loc+2],
        levels=[0],
        colors='k',
        linewidths=0.5,
        linestyles='--',
        zorder=0
    )

    # 绘制海床地形 (使用一维距离和深度)
    ax2.plot(plot_x, seafloor_depth, 'k-', linewidth=2, label='海床地形')
    ax2.fill_between(plot_x, seafloor_depth, max_y, color='lightgray', zorder=10)

    # 添加海平面线
    ax2.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

    # 设置坐标轴
    ax2.invert_yaxis()  # 深度向下增加
    ax2.set_ylim(max_y, 0)
    ax2.set_title('Summer', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.75)

    # 秋季
    ax3 = plt.subplot(223)
    # 创建正确的二维网格坐标
    X_grid, Y_grid = np.meshgrid(plot_x, autumn_data['y'][:max_y_loc+2], indexing='ij')  # 确保二维网格

    # 绘制水温数据 (使用二维坐标)
    mesh = ax3.contourf(
        X_grid, Y_grid, autumn_data['data'][:,:max_y_loc+2],  # 使用二维坐标和流速数据
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        levels=np.arange(vmin, vmax+0.05, 0.05),
        zorder=0
    )

    # 添加等值线
    ax3.contour(
        X_grid, Y_grid, autumn_data['data'][:,:max_y_loc+2],
        levels=[0],
        colors='k',
        linewidths=0.5,
        linestyles='--',
        zorder=0
    )

    # 绘制海床地形 (使用一维距离和深度)
    ax3.plot(plot_x, seafloor_depth, 'k-', linewidth=2, label='海床地形')
    ax3.fill_between(plot_x, seafloor_depth, max_y, color='lightgray', zorder=10)

    # 添加海平面线
    ax3.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

    # 设置坐标轴
    ax3.invert_yaxis()  # 深度向下增加
    ax3.set_ylim(max_y, 0)
    ax3.set_xlabel(x_label, fontsize=12)
    ax3.set_ylabel('Depth (m)', fontsize=12)
    ax3.set_title('Autumn', fontsize=12)
    ax3.grid(True, linestyle='--', alpha=0.75)

    # 冬季
    ax4 = plt.subplot(224)
    # 创建正确的二维网格坐标
    X_grid, Y_grid = np.meshgrid(plot_x, winter_data['y'][:max_y_loc+2], indexing='ij')  # 确保二维网格

    # 绘制水温数据 (使用二维坐标)
    mesh = ax4.contourf(
        X_grid, Y_grid, winter_data['data'][:,:max_y_loc+2],  # 使用二维坐标和流速数据
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        levels=np.arange(vmin, vmax+0.05, 0.05),
        zorder=0
    )

    # 添加等值线
    ax4.contour(
        X_grid, Y_grid, winter_data['data'][:,:max_y_loc+2],
        levels=[0],
        colors='k',
        linewidths=0.5,
        linestyles='--',
        zorder=0
    )

    # 绘制海床地形 (使用一维距离和深度)
    ax4.plot(plot_x, seafloor_depth, 'k-', linewidth=2, label='海床地形')
    ax4.fill_between(plot_x, seafloor_depth, max_y, color='lightgray', zorder=10)

    # 添加海平面线
    ax4.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

    # 设置坐标轴
    ax4.invert_yaxis()  # 深度向下增加
    ax4.set_ylim(max_y, 0)
    ax4.set_xlabel(x_label, fontsize=12)
    ax4.set_title('Winter', fontsize=12)
    ax4.grid(True, linestyle='--', alpha=0.75)

    # 在右侧添加一个共享的colorbar（垂直方向）
    cbar_ax = fig.add_axes([0.88, 0.08, 0.02, 0.82])  # [left, bottom, width, height]
    cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='vertical')
    cbar.set_label(r'Water nomrmal velocity $(m \cdot s^{-1})$', fontsize=12)

    plt.savefig('output/HYCOM/'+current_section+'/季节法向速度断面.png', dpi=600)
    plt.show()
    plt.close()