import glob
import netCDF4
from datetime import datetime,timedelta
import numpy as np
import matplotlib.pyplot as plt
import toml
from geopy.distance import geodesic as gd
import cmocean
plt.rcParams['font.family'] = 'Times New Roman, SimHei'
plt.rcParams['mathtext.fontset'] = 'stix'

# 计算断面法向量和断面的方位角
def calculate_normal_vector(lat1, lon1, lat2, lon2):  # 计算法向量与方位角
    # 计算中点纬度（用于投影）
    mid_lat = (lat1 + lat2) / 2
    # 计算经度差（考虑东西方向）
    dlon = lon2 - lon1
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360
    # 计算纬度差
    dlat = lat2 - lat1
    # 计算方位角
    azimuth = np.arctan2(dlon * np.cos(np.radians(mid_lat)), dlat)
    # 法向单位向量
    n_north = -np.sin(azimuth)
    n_east = np.cos(azimuth)
    # 归一化
    norm = np.sqrt(n_north ** 2 + n_east ** 2)
    # 法向单位向量（y前x后），起始到结束点的方位角
    if n_east < 0:
        return -n_north / norm, -n_east / norm, azimuth  # 将法向量改为东向
    else:
        return n_north / norm, n_east / norm, azimuth  # 本来是东向，不改

def calc_distance(x0,y0,x1,y1):
    return ((x1-x0)**2+(y1-y0)**2)**0.5

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

def interp_IDW(x, y, z, x0, y0):    # x,y自变量;z因变量;x0,y0插值点
    n_sta = len(x)
    dist = [] # 格点至所有站点的距离
    for s in range(n_sta):
        d = calc_distance(x[s], y[s], x0, y0)
        dist.append(d)
    wgt = 1.0 / np.power(dist, 2)
    wgt_sum = np.sum(wgt)
    arg_sum = np.sum(np.array(wgt) * np.array(z))
    result = arg_sum / wgt_sum
    return result

def get_distance(lat1, lon1, lat2, lon2):
    # 两点的经度和纬度
    p1 = (lat1, lon1)
    p2 = (lat2, lon2)

    d = gd(p1, p2).m	# 得到米为单位的结果
    return d

# 主程序
if __name__ == "__main__":
    water_temp_all_combined = []
    salinity_all_combined = []
    time_all_combined = []

    # 读取设置——注:用不到大区域故将lat_max等变量名给数据范围
    config = toml.load('config/config.toml')
    current_section = config['line_config']['current_section']['section']
    files = sorted(glob.glob(config['file_config']['HYCOM_reanalysis_ts']))
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

    # 创建存储数组
    water_temp_used = np.zeros((len(files),num_points,40),dtype=float)   # (time,points,levels)
    water_temp_used[:,:,:] = np.nan
    salinity_used = np.zeros((len(files),num_points,40),dtype=float)   # (time,points,levels)
    salinity_used[:,:,:] = np.nan

    # 设置绘图最深深度
    # max_y = (seafloor_depth.max()//100+1)*100 # 最大深度
    max_y = 1500

    # 计算每一时刻
    for i_file,file in enumerate(files):
        # === 2. 处理HYCOM流速数据 ===
        # 读取数据并筛选经纬度范围
        print("\n处理HYCOM流速数据...")
        hycom_data = netCDF4.Dataset(file, 'r')
        hycom_lon = hycom_data.variables['lon'][:]  # lon[-180. …… 180]
        hycom_lat = hycom_data.variables['lat'][:]  # lat[-80. …… 90]
        depth_levels = hycom_data.variables['depth'][:]
        print(f"深度层: {depth_levels}")
        hour = hycom_data.variables['time'][:]
        current_time = datetime(2000,1,1)+timedelta(hours=hour[0])
        current_time_str = current_time.strftime('%Y年%m月')
        time_all_combined.append(current_time_str)
        print('当前时间:', current_time_str)
        # (40,31,13) 层,纬度,经度
        water_temp   = hycom_data.variables['water_temp'][0, :, :, :]   # (time,depth,lat,lon)
        salinity   = hycom_data.variables['salinity'][0, :, :, :]   # (time,depth,lat,lon)

        # 获取对应点的速度数据
        for i_point in range(len(points)):
            for i_lat in range(len(hycom_lat)):
                if hycom_lat[i_lat]<points[i_point,0]<hycom_lat[i_lat+1]:
                    lat_loc = i_lat
                    break
            for i_lon in range(len(hycom_lon)):
                if hycom_lon[i_lon]<points[i_point,1]<hycom_lon[i_lon+1]:
                    lon_loc = i_lon
                    break

            distances_used = []
            lat_used = []
            lon_used = []
            for i_lat in range(lat_loc-2,lat_loc+2):
                for i_lon in range(lon_loc-2,lon_loc+2):
                    lat_used.append(i_lat)
                    lon_used.append(i_lon)
                    distances_used.append(get_distance(hycom_lat[i_lat],hycom_lon[i_lon],points[i_point][0],points[i_point][1]))
            min_distance = np.argmin(distances_used)

            water_temp_used[i_file,i_point,:] = water_temp[:,lat_used[min_distance],lon_used[min_distance]]
            salinity_used[i_file,i_point,:] = salinity[:,lat_used[min_distance],lon_used[min_distance]]

        # === 3. 存储新的坐标、速度数据 ===
        print("\n存储完整的绘图变量...")
        # 存储深度
        plot_depth = np.zeros((num_points,len(depth_levels)),dtype=float)
        for i_point in range(num_points):
            plot_depth[i_point,:] = depth_levels[:]
        # 存储流速
        water_temp_used[water_temp_used == -30000] = np.nan
        salinity_used[salinity_used == -30000] = np.nan

        # === 4. 计算水平距离 ===
        print("\n计算水平距离...")
        # 计算点间距离（km）
        segment_distances = np.zeros(num_points - 1)
        for i in range(num_points - 1):
            p1 = (lat_list[i], lon_list[i])
            p2 = (lat_list[i + 1], lon_list[i + 1])
            segment_distances[i] = gd(p1, p2).km
        # 转换为米
        point_widths_m = segment_distances * 1000.0
        plot_grid = []
        for i_point in range(num_points):
            plot_grid.append(sum(point_widths_m[:i_point]))
        plot_grid = np.array(plot_grid)

        # === 5. 得到温盐数据 ===
        print("\n得到温盐数据...")
        valid_water_temp = np.zeros((num_points, len(depth_levels)))
        valid_salinity = np.zeros((num_points, len(depth_levels)))
        valid_water_temp[:, :] = np.nan
        valid_salinity[:, :] = np.nan
        for i_point in range(num_points):
            i_depth = 0
            i_bottom = np.sum(depth_levels < seafloor_depth[i_point])
            for i_depth in range(depth_levels.shape[0]):
                if i_depth <= i_bottom-1:
                    valid_salinity[i_point, i_depth] = salinity_used[i_file,i_point, i_depth]
                    valid_water_temp[i_point, i_depth] = water_temp_used[i_file,i_point,i_depth]
                else:
                    valid_salinity[i_point, i_depth] = np.nan
                    valid_water_temp[i_point, i_depth] = np.nan

        water_temp_all_combined.append(valid_water_temp)
        salinity_all_combined.append(valid_salinity)

        # === 8. 创建网格用于绘图 ===
        print("\n创建绘图网格...")

        # 计算累积距离
        distances = np.zeros(num_points)
        for i in range(1, num_points):
            distances[i] = distances[i - 1] + segment_distances[i - 1]

        # 深度网格
        X_grid, Y_grid = np.meshgrid(distances, depth_levels, indexing='ij')

        # === 9. 绘制截面水温图 ===
        print("\n绘制截面水温图...")
        plt.figure(figsize=(14, 8))

        # 创建正确的二维网格坐标
        X_grid, Y_grid = np.meshgrid(plot_grid, depth_levels, indexing='ij')  # 确保二维网格

        # 创建颜色映射
        cmap = cmocean.cm.balance
        cmap.set_bad('lightgray', 1.0)

        # 计算颜色范围 (基于截面水温)
        vmax = np.nanmax(valid_water_temp)
        vmin = np.nanmin(valid_water_temp)

        # 绘制水温数据 (使用二维坐标)
        mesh = plt.contourf(
            X_grid/1000., Y_grid, valid_water_temp,  # 使用二维坐标和水温数据
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=50,  # 明确指定层级数
            zorder=0
        )

        # 添加零线
        plt.contour(
            X_grid/1000., Y_grid, valid_water_temp,
            levels=[0],
            colors='k',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 绘制海床地形 (使用一维距离和深度)
        plt.plot(plot_grid/1000., seafloor_depth, 'k-', linewidth=2, label='海床地形')
        plt.fill_between(plot_grid/1000., seafloor_depth, max_y, color='lightgray', zorder=10)

        # 添加海平面线
        plt.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

        # 添加颜色条
        cbar = plt.colorbar(mesh, label='水温 (℃)')
        cbar.set_label('水温 (℃)', fontsize=12)

        # 设置坐标轴
        plt.gca().invert_yaxis()  # 深度向下增加
        plt.ylim(max_y, 0)
        plt.xlabel('沿断面距离 (km)', fontsize=12)
        plt.ylabel('深度 (m)', fontsize=12)
        plt.title(current_time_str+'水温截面图', fontsize=14)

        plt.grid(True, linestyle='--', alpha=0.75)
        plt.tight_layout()
        plt.savefig('output/HYCOM/'+current_section+'/'+current_time_str+'water_temp_section.png', dpi=600)
        plt.show()
        plt.close()

        # === 10. 绘制截面盐度图 ===
        print("\n绘制截面盐度图...")
        plt.figure(figsize=(14, 8))

        # 创建正确的二维网格坐标
        X_grid, Y_grid = np.meshgrid(plot_grid, depth_levels, indexing='ij')  # 确保二维网格

        # 创建颜色映射
        cmap = cmocean.cm.balance
        cmap.set_bad('lightgray', 1.0)

        # 计算颜色范围 (基于法向流速)
        vmax = np.nanmax(valid_salinity)
        vmin = np.nanmin(valid_salinity)

        # 绘制水温数据 (使用二维坐标)
        mesh = plt.contourf(
            X_grid/1000., Y_grid, valid_salinity,  # 使用二维坐标和流速数据
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=50,  # 明确指定层级数
            zorder=0
        )

        # 添加零线
        plt.contour(
            X_grid/1000., Y_grid, valid_salinity,
            levels=[0],
            colors='k',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 绘制海床地形 (使用一维距离和深度)
        plt.plot(plot_grid/1000., seafloor_depth, 'k-', linewidth=2, label='海床地形')
        plt.fill_between(plot_grid/1000., seafloor_depth, max_y, color='lightgray', zorder=10)

        # 添加海平面线
        plt.axhline(0, color='blue', linestyle='-', linewidth=1, alpha=0.7, zorder=1)

        # 添加颜色条
        cbar = plt.colorbar(mesh, label='盐度 (PSU)')
        cbar.set_label('盐度 (PSU)', fontsize=12)

        # 设置坐标轴
        plt.gca().invert_yaxis()  # 深度向下增加
        plt.ylim(max_y, 0)
        plt.xlabel('沿断面距离 (km)', fontsize=12)
        plt.ylabel('深度 (m)', fontsize=12)
        plt.title(current_time_str+'盐度截面图', fontsize=14)

        plt.grid(True, linestyle='--', alpha=0.75)
        plt.tight_layout()
        plt.savefig('output/HYCOM/'+current_section+'/'+current_time_str+'salinity_section.png', dpi=600)
        plt.show()
        plt.close()

    # 输出数据
    water_temp_all_combined = np.array(water_temp_all_combined)
    salinity_all_combined = np.array(salinity_all_combined)
    spring_water_temp = []
    summer_water_temp = []
    autumn_water_temp = []
    winter_water_temp = []
    for i in range(water_temp_all_combined.shape[0]):
        if '12月' in time_all_combined[i] or '01月' in time_all_combined[i] or '02月' in time_all_combined[i]:   # 冬天
            winter_water_temp.append(water_temp_all_combined[i])
        elif '3月' in time_all_combined[i] or '4月' in time_all_combined[i] or '5月' in time_all_combined[i]:  # 春天
            spring_water_temp.append(water_temp_all_combined[i])
        elif '6月' in time_all_combined[i] or '7月' in time_all_combined[i] or '8月' in time_all_combined[i]:  # 夏天
            summer_water_temp.append(water_temp_all_combined[i])
        elif '9月' in time_all_combined[i] or '10月' in time_all_combined[i] or '11月' in time_all_combined[i]:# 秋天
            autumn_water_temp.append(water_temp_all_combined[i])
    spring_water_temp = np.array(spring_water_temp)
    summer_water_temp = np.array(summer_water_temp)
    autumn_water_temp = np.array(autumn_water_temp)
    winter_water_temp = np.array(winter_water_temp)
    np.savez('output/HYCOM/'+current_section+'/spring_water_temp_section',x=distances, y=depth_levels, t=time_all_combined, data=spring_water_temp)
    np.savez('output/HYCOM/'+current_section+'/summer_water_temp_section',x=distances, y=depth_levels, t=time_all_combined, data=summer_water_temp)
    np.savez('output/HYCOM/'+current_section+'/autumn_water_temp_section',x=distances, y=depth_levels, t=time_all_combined, data=autumn_water_temp)
    np.savez('output/HYCOM/'+current_section+'/winter_water_temp_section',x=distances, y=depth_levels, t=time_all_combined, data=winter_water_temp)
    ave_spring_water_temp = np.nanmean(spring_water_temp[:],axis=0)
    ave_summer_water_temp = np.nanmean(summer_water_temp[:],axis=0)
    ave_autumn_water_temp = np.nanmean(autumn_water_temp[:],axis=0)
    ave_winter_water_temp = np.nanmean(winter_water_temp[:],axis=0)
    np.savez('output/HYCOM/'+current_section+'/ave_spring_water_temp_section',x=distances, y=depth_levels, data=ave_spring_water_temp)
    np.savez('output/HYCOM/'+current_section+'/ave_summer_water_temp_section',x=distances, y=depth_levels, data=ave_summer_water_temp)
    np.savez('output/HYCOM/'+current_section+'/ave_autumn_water_temp_section',x=distances, y=depth_levels, data=ave_autumn_water_temp)
    np.savez('output/HYCOM/'+current_section+'/ave_winter_water_temp_section',x=distances, y=depth_levels, data=ave_winter_water_temp)
    
    winter_salinity = []
    spring_salinity = []
    summer_salinity = []
    autumn_salinity = []
    for i in range(salinity_all_combined.shape[0]):
        if '12月' in time_all_combined[i] or '1月' in time_all_combined[i] or '2月' in time_all_combined[i]:   # 冬天
            winter_salinity.append(salinity_all_combined[i])
        elif '3月' in time_all_combined[i] or '4月' in time_all_combined[i] or '5月' in time_all_combined[i]:  # 春天
            spring_salinity.append(salinity_all_combined[i])
        elif '6月' in time_all_combined[i] or '7月' in time_all_combined[i] or '8月' in time_all_combined[i]:  # 夏天
            summer_salinity.append(salinity_all_combined[i])
        elif '9月' in time_all_combined[i] or '10月' in time_all_combined[i] or '11月' in time_all_combined[i]:# 秋天
            autumn_salinity.append(salinity_all_combined[i])
    spring_salinity = np.array(spring_salinity)
    summer_salinity = np.array(summer_salinity)
    autumn_salinity = np.array(autumn_salinity)
    winter_salinity = np.array(winter_salinity)
    np.savez('output/HYCOM/'+current_section+'/spring_salinity_section',x=distances, y=depth_levels, t=time_all_combined, data=spring_salinity)
    np.savez('output/HYCOM/'+current_section+'/summer_salinity_section',x=distances, y=depth_levels, t=time_all_combined, data=summer_salinity)
    np.savez('output/HYCOM/'+current_section+'/autumn_salinity_section',x=distances, y=depth_levels, t=time_all_combined, data=autumn_salinity)
    np.savez('output/HYCOM/'+current_section+'/winter_salinity_section',x=distances, y=depth_levels, t=time_all_combined, data=winter_salinity)
    ave_spring_salinity = np.nanmean(spring_salinity[:],axis=0)
    ave_summer_salinity = np.nanmean(summer_salinity[:],axis=0)
    ave_autumn_salinity = np.nanmean(autumn_salinity[:],axis=0)
    ave_winter_salinity = np.nanmean(winter_salinity[:],axis=0)
    np.savez('output/HYCOM/'+current_section+'/ave_spring_salinity_section',x=distances, y=depth_levels, data=ave_spring_salinity)
    np.savez('output/HYCOM/'+current_section+'/ave_summer_salinity_section',x=distances, y=depth_levels, data=ave_summer_salinity)
    np.savez('output/HYCOM/'+current_section+'/ave_autumn_salinity_section',x=distances, y=depth_levels, data=ave_autumn_salinity)
    np.savez('output/HYCOM/'+current_section+'/ave_winter_salinity_section',x=distances, y=depth_levels, data=ave_winter_salinity)