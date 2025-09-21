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
    speed_all_combined = []
    flux_all_combined = []
    time_all_combined = []

    # 读取设置——注:用不到大区域故将lat_max等变量名给数据范围
    config = toml.load('config/config.toml')
    current_section = config['line_config']['current_section']['section']
    files = sorted(glob.glob(config['file_config']['HYCOM_reanalysis_uv']))
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
    u_used = np.zeros((len(files),num_points,40),dtype=float)   # (time,points,levels)
    u_used[:,:,:] = np.nan
    v_used = np.zeros((len(files),num_points,40),dtype=float)   # (time,points,levels)
    v_used[:,:,:] = np.nan

    # 设置绘图最深深度
    max_y = (seafloor_depth.max()//100+1)*100 # 最大深度

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
        print('当前时间:', current_time_str)
        # (40,31,13) 层,纬度,经度
        u_vel   = hycom_data.variables['water_u'][0, :, :, :]   # (time,depth,lat,lon)
        v_vel   = hycom_data.variables['water_v'][0, :, :, :]   # (time,depth,lat,lon)

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

            u_used[i_file,i_point,:] = u_vel[:,lat_used[min_distance],lon_used[min_distance]]
            v_used[i_file,i_point,:] = u_vel[:,lat_used[min_distance],lon_used[min_distance]]

        # === 3. 存储新的坐标、速度数据 ===
        print("\n存储完整的绘图变量...")
        # 存储深度
        plot_depth = np.zeros((num_points,len(depth_levels)),dtype=float)
        for i_point in range(num_points):
            plot_depth[i_point,:] = depth_levels[:]
        # 存储流速
        u_used[u_used == 30000] = np.nan
        v_used[v_used == 30000] = np.nan
        plot_u = u_used[i_file,:,:]
        plot_v = v_used[i_file,:,:]

        # === 4. 计算法向流速 ===
        print("\n计算法向流速分量...")
        # 计算断面法向量和方位角
        north_comp, east_comp, azimuth = calculate_normal_vector(point_start_loc[0], point_start_loc[1],point_end_loc[0], point_end_loc[1])
        print(f"断面法向量 - 北分量: {north_comp:.6f}, 东分量: {east_comp:.6f}")
        print(f"断面方位角: {np.degrees(azimuth):.2f}°")
        # 计算径向单位向量（沿断面方向）
        r_north = np.cos(azimuth)
        r_east = np.sin(azimuth)
        print(f"断面径向向量 - 北分量: {r_north:.6f}, 东分量: {r_east:.6f}")
        # 计算流速分量
        normal_velocities = np.zeros(plot_u.shape,dtype=float)  # 垂直于断面方向
        # 计算各层流速分量
        for i in range(num_points):
            for k in range(len(depth_levels)):
                # 获取流速分量
                u = plot_u[i, k]  # 东向分量
                v = plot_v[i, k]  # 北向分量
                # 计算法向分量: V_normal = u * n_east + v * n_north
                normal_velocities[i, k] = u * east_comp + v * north_comp
        # 删除异常流速值
        normal_velocities[abs(normal_velocities)>10.0] = np.nan
        # 打印原始流速范围
        print(f"法向流速范围: {np.nanmin(normal_velocities):.6f} 到 {np.nanmax(normal_velocities):.6f} m/s")

        # === 5. 计算水平距离 ===
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
        actual_distance = []
        for i_point in range(num_points):
            if i_point == 0:
                actual_distance.append(point_widths_m[0]/2)
            elif i_point == num_points-1:
                actual_distance.append(point_widths_m[-1]/2)
            else:
                actual_distance.append(point_widths_m[i_point-1]/2+point_widths_m[i_point]/2)
        actual_distance = np.array(actual_distance)

        # === 6. 计算面积和通量(m²,m³/s) ===
        print("\n计算面积和通量...")
        valid_areas = np.zeros((num_points, len(depth_levels)))
        valid_speed = np.zeros((num_points, len(depth_levels)))
        valid_flux  = np.zeros((num_points, len(depth_levels)))
        valid_areas[:, :] = np.nan
        valid_speed[:, :] = np.nan
        valid_flux[:, :]  = np.nan
        for i_point in range(num_points):
            i_depth = 0
            i_bottom = np.sum(depth_levels < seafloor_depth[i_point])
            for i_depth in range(depth_levels.shape[0]):
                if i_depth == 0:
                    valid_areas[i_point, i_depth] = actual_distance[i_point]*(depth_levels[i_depth+1]-depth_levels[i_depth])/2  # 下层的一半
                    valid_speed[i_point, i_depth] = normal_velocities[i_point, i_depth]
                elif i_depth < i_bottom-1:
                    valid_areas[i_point, i_depth] = actual_distance[i_point]*((depth_levels[i_depth+1]-depth_levels[i_depth])/2+(depth_levels[i_depth]-depth_levels[i_depth-1])/2)  # 下层的一半+上层的一半
                    valid_speed[i_point, i_depth] = normal_velocities[i_point, i_depth]
                elif i_depth == i_bottom-1:
                    valid_areas[i_point, i_depth] = depth_levels[i_depth]-depth_levels[i_depth-1]+actual_distance[i_point]*(seafloor_depth[i_point]-depth_levels[i_depth])/2    # 上层的一半+下层的全部
                    valid_speed[i_point, i_depth] = normal_velocities[i_point, i_depth]
                else:
                    valid_areas[i_point, i_depth] = np.nan
                    valid_speed[i_point, i_depth] = np.nan

        speed_all_combined.append(valid_speed)

        # 验证有效面积
        max_area = np.nanmax(valid_areas)
        min_area = np.nanmin(valid_areas[valid_areas > 0])
        print(f"最大有效面积: {max_area:.2f} m²")
        print(f"最小有效面积: {min_area:.2f} m²")

        # === 7. 计算海水通量 ===
        print("\n计算海水通量...")

        # 计算通量 (Sv)
        current_flux = valid_speed * valid_areas
        total_flux_sv = np.nansum(current_flux)/ 1e6
        flux_all_combined.append(total_flux_sv)

        print(f"断面海水总通量: {total_flux_sv:.6f} Sv")

        # === 8. 创建网格用于绘图 ===
        print("\n创建绘图网格...")

        # 计算累积距离
        distances = np.zeros(num_points)
        for i in range(1, num_points):
            distances[i] = distances[i - 1] + segment_distances[i - 1]

        # 深度网格
        X_grid, Y_grid = np.meshgrid(distances, depth_levels, indexing='ij')

        # === 9. 绘制法向流速分量图 ===
        print("\n绘制法向流速分量图...")
        plt.figure(figsize=(14, 8))

        # 创建正确的二维网格坐标
        X_grid, Y_grid = np.meshgrid(plot_grid, depth_levels, indexing='ij')  # 确保二维网格

        # 创建颜色映射
        cmap = cmocean.cm.balance
        cmap.set_bad('lightgray', 1.0)

        # 计算颜色范围 (基于法向流速)
        vmax = np.nanmax(abs(normal_velocities))
        vmin = -vmax

        # 绘制法向流速数据 (使用二维坐标)
        mesh = plt.contourf(
            X_grid/1000., Y_grid, normal_velocities,  # 使用二维坐标和流速数据
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=50,  # 明确指定层级数
            zorder=0
        )

        # 添加零线
        plt.contour(
            X_grid/1000., Y_grid, normal_velocities,
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
        cbar = plt.colorbar(mesh, label='法向流速 (m/s)')
        cbar.set_label('法向流速 (m/s)', fontsize=12)

        # 设置坐标轴
        plt.gca().invert_yaxis()  # 深度向下增加
        plt.ylim(max_y, 0)
        plt.xlabel('沿断面距离 (km)', fontsize=12)
        plt.ylabel('深度 (m)', fontsize=12)
        plt.title(current_time_str+'法向流速分量截面图', fontsize=14)

        # 添加图例
        plt.legend(loc='lower right')

        plt.grid(True, linestyle='--', alpha=0.75)
        plt.tight_layout()
        plt.savefig('output/HYCOM/'+current_section+'/'+current_time_str+'normal_velocity_section.png', dpi=300)
        plt.show()

        # 存储通量数据入数组
        time_all_combined.append(current_time_str)

    # 输出数据
    speed_all_combined = np.array(speed_all_combined)
    flux_all_combined = np.array(flux_all_combined)
    # 保存通量
    np.savez('output/HYCOM/'+current_section+'/volume_transport',t=time_all_combined,data=flux_all_combined)
    # 保存流速
    spring_speed = []
    summer_speed = []
    autumn_speed = []
    winter_speed = []
    for i in range(speed_all_combined.shape[0]):
        if '12月' in time_all_combined[i] or '01月' in time_all_combined[i] or '02月' in time_all_combined[i]:   # 冬天
            winter_speed.append(speed_all_combined[i])
        elif '3月' in time_all_combined[i] or '4月' in time_all_combined[i] or '5月' in time_all_combined[i]:  # 春天
            spring_speed.append(speed_all_combined[i])
        elif '6月' in time_all_combined[i] or '7月' in time_all_combined[i] or '8月' in time_all_combined[i]:  # 夏天
            summer_speed.append(speed_all_combined[i])
        elif '9月' in time_all_combined[i] or '10月' in time_all_combined[i] or '11月' in time_all_combined[i]:# 秋天
            autumn_speed.append(speed_all_combined[i])
    spring_speed = np.array(spring_speed)
    summer_speed = np.array(summer_speed)
    autumn_speed = np.array(autumn_speed)
    winter_speed = np.array(winter_speed)
    np.savez('output/HYCOM/'+current_section+'/spring_velocity_section',x=distances, y=depth_levels, t=time_all_combined, data=spring_speed)
    np.savez('output/HYCOM/'+current_section+'/summer_velocity_section',x=distances, y=depth_levels, t=time_all_combined, data=summer_speed)
    np.savez('output/HYCOM/'+current_section+'/autumn_velocity_section',x=distances, y=depth_levels, t=time_all_combined, data=autumn_speed)
    np.savez('output/HYCOM/'+current_section+'/winter_velocity_section',x=distances, y=depth_levels, t=time_all_combined, data=winter_speed)
    ave_spring_speed = np.nanmean(spring_speed[:],axis=0)
    ave_summer_speed = np.nanmean(summer_speed[:],axis=0)
    ave_autumn_speed = np.nanmean(autumn_speed[:],axis=0)
    ave_winter_speed = np.nanmean(winter_speed[:],axis=0)
    np.savez('output/HYCOM/'+current_section+'/ave_spring_velocity_section',x=distances, y=depth_levels, data=ave_spring_speed)
    np.savez('output/HYCOM/'+current_section+'/ave_summer_velocity_section',x=distances, y=depth_levels, data=ave_summer_speed)
    np.savez('output/HYCOM/'+current_section+'/ave_autumn_velocity_section',x=distances, y=depth_levels, data=ave_autumn_speed)
    np.savez('output/HYCOM/'+current_section+'/ave_winter_velocity_section',x=distances, y=depth_levels, data=ave_winter_speed)