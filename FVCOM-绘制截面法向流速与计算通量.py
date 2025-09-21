import glob
import netCDF4
from datetime import datetime,timedelta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
import toml
from scipy.interpolate import RegularGridInterpolator
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

# 主程序
if __name__ == "__main__":
    flux_all_combined = []
    time_all_combined = []
    # files = glob.glob(r'D:\Documents\code\Pycharm\FVCOM_analysis\daily_average\daily_average\*.nc')
    file = pathlib.Path(r'D:\Documents\code\Pycharm\FVCOM_analysis\kuroshio_2024_0717\kuroshio_2024_0717.nc')

    # 预处理
    max_y = 1200 # 最大深度

    # 读取设置——注:用不到大区域故将lat_max等变量名给数据范围
    config = toml.load('config/config.toml')
    line_start_lat = config['line_config']['start_lat']
    line_start_lon = config['line_config']['start_lon']
    line_end_lat = config['line_config']['end_lat']
    line_end_lon = config['line_config']['end_lon']
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

    # 计算每一时刻
    # === 2. 处理fvcom流速数据 ===
    # 读取数据并筛选经纬度范围
    print("\n处理fvcom流速数据...")
    fvcom_data = netCDF4.Dataset(file, 'r')
    fvcom_lonc = fvcom_data.variables['lonc'][:]        # (441027)      网格
    fvcom_latc = fvcom_data.variables['latc'][:]        # (441027)      网格
    fvcom_siglay = fvcom_data.variables['siglay'][:]    # (68,441027)
    fvcom_u = fvcom_data.variables['u'][:]              # (24,68,441027)
    fvcom_v = fvcom_data.variables['v'][:]              # (24,68,441027)
    fvcom_times = fvcom_data.variables['Times'][:]      # (24,26)
    fvcom_lat = fvcom_data.variables['lat'][:]          # (226031)      仅用于水深插值
    fvcom_lon = fvcom_data.variables['lon'][:]          # (226031)      仅用于水深插值
    fvcom_h = fvcom_data.variables['h'][:]              # (226031)      仅用于水深插值
    fvcom_zeta = fvcom_data.variables['zeta'][:]        # (24,226031)   仅用于水深插值
    times_all = []
    for i_time in range(len(fvcom_times)):
        temp = ''
        for j in range(len(fvcom_times[0])):
            temp+=fvcom_times[i_time,j].decode('utf-8')
        times_all.append((pd.to_datetime(temp)+timedelta(hours=8)).strftime('%Y-%m-%d %H时'))
    cell_used = []
    # 筛选区域内及附近网格
    for i_cell in range(len(fvcom_lonc)):
        if lon_min-1 < fvcom_lonc[i_cell] < lon_max+1 and lat_min-1 < fvcom_latc[i_cell] < lat_max+1:
            cell_used.append(i_cell)
    fvcom_data.close()
    # 筛选计算绘图要用的 latc,lonc,siglay
    latc_used = fvcom_latc[cell_used]
    lonc_used = fvcom_lonc[cell_used]
    siglay_used = fvcom_siglay[:,cell_used]
    depth = np.zeros_like(fvcom_zeta,dtype=float)
    depth[:,:] = np.nan
    for time_idx in range(len(fvcom_times)):
        depth[time_idx,:] = fvcom_zeta[time_idx,:]+fvcom_h[:]

    # 逐小时处理
    for time_idx in range(len(fvcom_times)):
        # 筛选计算绘图要用的 u, v, depth
        u_used = fvcom_u[time_idx,:,cell_used]
        v_used = fvcom_v[time_idx,:,cell_used]
        depth_temp = depth[time_idx,:]      # 注 还需要插值到网格
        # 创建水深插值函数并插值
        depth_interp = RegularGridInterpolator(
            (fvcom_lat, fvcom_lon),
            depth_temp,
            method='linear',
            bounds_error=False,
            fill_value=np.nan
        )
        depth_used = depth_interp(points)

        # 创建插值网格
        u_grid = np.zeros((num_points, len(depth_levels)))
        v_grid = np.zeros((num_points, len(depth_levels)))

        #每个深度插值
        for depth_idx in range(len(depth_levels)):
            # 创建速度插值函数并插值
            u_interp = RegularGridInterpolator(
                (lat_used, lon_used),
                u_used[depth_idx],
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            v_interp = RegularGridInterpolator(
                (lat_used, lon_used),
                v_used[depth_idx],
                method='linear',
                bounds_error=False,
                fill_value=np.nan
            )
            u_grid[:,depth_idx] = u_interp(points)
            v_grid[:,depth_idx] = v_interp(points)

        # 检查深度与流速匹配关系
        for point_idx in range(num_points):  # 点
            depth_current_point = seafloor_depth[point_idx]
            # 当当前深度大于该点深度时，全部为np.nan
            depth_current_idx = np.argwhere(depth_levels > depth_current_point)[0][0]
            u_grid[point_idx, depth_current_idx:] = np.nan
            v_grid[point_idx, depth_current_idx:] = np.nan
            # 判断底部流速是否符合常理
            for depth_idx in range(len(depth_levels)):
                if np.isnan(u_grid[point_idx, depth_idx]):
                    continue
                elif np.abs(u_grid[point_idx, depth_idx])>2 or np.abs(v_grid[point_idx, depth_idx])>2:  # 插值处理异常数据
                    try:
                        u_grid[point_idx, depth_idx] = np.interp(depth_levels[depth_idx], depth_levels[:depth_idx], u_grid[point_idx, :depth_idx])
                    except ValueError:
                        u_grid[point_idx, depth_idx] = -99999.
                    try:
                        v_grid[point_idx, depth_idx] = np.interp(depth_levels[depth_idx], depth_levels[:depth_idx], v_grid[point_idx, :depth_idx])
                    except ValueError:
                        v_grid[point_idx, depth_idx] = -99999.
                else:
                    continue
        # 处理整层异常数据，从横向插值
        for depth_idx in range(len(depth_levels)):
            err_surface_idx = []
            cor_surface_idx = []
            cor_u_spd = []
            cor_v_spd = []
            for point_idx in range(num_points):
                if abs(u_grid[point_idx, depth_idx])>2 or abs(v_grid[point_idx, depth_idx])>2:    # 异常数据
                    err_surface_idx.append(point_idx)
                else:                                                   # 正常数据
                    cor_surface_idx.append(point_idx)
                    cor_u_spd.append(u_grid[point_idx, depth_idx])
                    cor_v_spd.append(v_grid[point_idx, depth_idx])
            err_surface_idx = np.array(err_surface_idx)
            cor_surface_idx = np.array(cor_surface_idx)
            cor_u_spd = np.array(cor_u_spd)
            cor_v_spd = np.array(cor_v_spd)
            corrected_surface_u = np.interp(err_surface_idx, cor_surface_idx, cor_u_spd)
            corrected_surface_v = np.interp(err_surface_idx, cor_surface_idx, cor_v_spd)
            for i_idx,idx in enumerate(err_surface_idx):
                u_grid[idx, depth_idx] = corrected_surface_u[i_idx]
                v_grid[idx, depth_idx] = corrected_surface_v[i_idx]

        # === 3. 存储新的坐标、速度数据 ===
        print("\n存储完整的绘图变量...")
        # 存储深度
        plot_depth = np.zeros((num_points,len(depth_levels)),dtype=float)
        plot_depth[:,:] = np.nan
        plot_seafloor_mark = np.zeros(num_points,dtype=int)
        for i_point in range(num_points):
            for i_depth in range(len(depth_levels)):
                # 判断深度是否大于海床
                if i_depth==0:
                    plot_depth[i_point, i_depth] = depth_levels[i_depth]
                elif depth_levels[i_depth] < seafloor_depth[i_point]:
                    plot_depth[i_point, i_depth] = depth_levels[i_depth]    # 小于则记录实际深度
                elif depth_levels[i_depth]>= seafloor_depth[i_point] and depth_levels[i_depth-1] < seafloor_depth[i_point]:
                    plot_depth[i_point,i_depth] = seafloor_depth[i_point]   # 大于则记录海床深度
                    plot_seafloor_mark[i_point] = i_depth                   # 记录每个点的海床深度下标
                else:
                    continue
        # 存储流速
        plot_u = np.zeros(plot_depth.shape,dtype=float)
        plot_v = np.zeros(plot_depth.shape,dtype=float)
        plot_u[:,:] = np.nan
        plot_v[:,:] = np.nan
        for i_point in range(num_points):
            for i_depth in range(len(depth_levels)):
                if plot_seafloor_mark[i_point] > i_depth:       # 没到深度则记录流速
                    plot_u[i_point, i_depth] = u_grid[i_point, i_depth]
                    plot_v[i_point, i_depth] = v_grid[i_point, i_depth]
                elif plot_seafloor_mark[i_point] == i_depth:    # 达到深度则记录底部流速
                    plot_u[i_point, i_depth] = u_bottom_grid[i_point]
                    plot_v[i_point, i_depth] = v_bottom_grid[i_point]
                else:                                           # 超过深度则缺测
                    plot_u[i_point, i_depth] = np.nan
                    plot_v[i_point, i_depth] = np.nan
        # 处理部分缺测情况
        for i_point in range(num_points):
            max_nonnan = np.nan
            err_idx = []
            cor_idx = []
            cor_depth = []
            cor_u_spd = []
            cor_v_spd = []
            for i_depth in range(plot_depth.shape[1]-1,-1,-1):                                                  # range(0,40,1)->range(39,-1,-1)
                if (not np.isnan(plot_u[i_point, i_depth])) or (not np.isnan(plot_v[i_point, i_depth])):        # 最大不缺测的记录位置
                    max_nonnan = i_depth
                    break
            for i_depth in range(max_nonnan):
                if np.isnan(plot_v[i_point, i_depth]):                                                          # 记录海床深度以上缺测的下标
                    err_idx.append(i_depth)
                else:
                    cor_idx.append(i_depth)
                    cor_depth.append(plot_depth[i_point, i_depth])
                    cor_u_spd.append(plot_u[i_point, i_depth])
                    cor_v_spd.append(plot_v[i_point, i_depth])
            err_idx = np.array(err_idx)
            if len(err_idx) > 0:
                temp1 = np.interp(err_idx, cor_depth, cor_u_spd)
                temp2 = np.interp(err_idx, cor_depth, cor_v_spd)
                for i,idx in enumerate(err_idx):
                    plot_u[i_point,idx] = temp1[i]
                    plot_v[i_point,idx] = temp2[i]

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
        normal_velocities = np.zeros(plot_u.shape)  # 垂直于断面方向
        # 计算各层流速分量
        for i in range(num_points):
            for k in range(len(depth_levels)):
                # 获取流速分量
                u = plot_u[i, k]  # 东向分量
                v = plot_v[i, k]  # 北向分量
                # 计算法向分量: V_normal = u * n_east + v * n_north
                normal_velocities[i, k] = u * east_comp + v * north_comp
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

        # === 6. 计算面积和通量(m²,m³/s) ===
        print("\n计算面积和通量...")
        valid_areas = np.zeros((num_points-1, len(depth_levels)-1))
        valid_speed = np.zeros((num_points-1, len(depth_levels)-1))
        valid_flux  = np.zeros((num_points-1, len(depth_levels)-1))
        valid_areas[:, :] = np.nan
        valid_speed[:, :] = np.nan
        valid_flux[:, :]  = np.nan
        for i_point in range(num_points-1):
            i_depth = 0
            while i_depth < len(depth_levels)-1:
                top_left = plot_depth[i_point, i_depth]
                top_right = plot_depth[i_point+1, i_depth]
                bottom_left = plot_depth[i_point, i_depth+1]
                bottom_right = plot_depth[i_point+1, i_depth+1]
                # 判断是否为底部
                if np.isnan(plot_depth[i_point, i_depth+2]) or np.isnan(plot_depth[i_point+1, i_depth+2]):  # 是底部
                    i_depth_temp = i_depth
                    if np.isnan(plot_depth[i_point, i_depth_temp+2]):                                       # 若左侧底部
                        while not np.isnan(plot_depth[i_point+1, i_depth_temp+2]):                          # 找到右侧底部
                            i_depth_temp += 1
                        # i_depth为左侧深度,i_depth_temp为右侧深度(右侧更深)
                        bottom_right = plot_depth[i_point+1, i_depth_temp]
                        valid_areas[i_point, i_depth] = ((bottom_left - top_left)+(bottom_right - top_right))*point_widths_m[i_point]/2         # m²
                        # 计算形心位置
                        a = bottom_left - top_left                                                          # 上底
                        b = bottom_right - top_right                                                        # 下底
                        h = point_widths_m[i_point]                                                         # 高
                        x = []
                        y = []
                        z = []
                        for j_depth in range(i_depth,i_depth+2):                                            # 左侧,上底
                            x.append(plot_depth[i_point,j_depth]-plot_depth[i_point,i_depth])
                            y.append(plot_grid[i_point+1]-plot_grid[i_point])
                            z.append(normal_velocities[i_point,j_depth])
                        for j_depth in range(i_depth,i_depth_temp+1):                                       # 右侧,下底
                            x.append(plot_depth[i_point+1,j_depth]-plot_depth[i_point+1,i_depth])
                            y.append(0)
                            z.append(normal_velocities[i_point+1,j_depth])
                        x0 = (2*a+b)*h/3/(a+b)
                        y0 = (a**2+b**2+a*b)/3/(a+b)
                        # 插值到形心位置
                        valid_speed[i_point, i_depth] = interp_IDW(x,y,z,x0,y0)
                    else:                                                                                   # 右侧底部
                        while not np.isnan(plot_depth[i_point, i_depth_temp+2]):                            # 找到左侧底部
                            i_depth_temp += 1
                        # i_depth为右侧深度,i_depth_temp为左侧深度(左侧更深)
                        bottom_right = plot_depth[i_point, i_depth_temp]
                        # 计算形心位置
                        a = bottom_right - top_right                                                        # 上底
                        b = bottom_left - top_left                                                          # 下底
                        h = point_widths_m[i_point]                                                         # 高
                        x = []
                        y = []
                        z = []
                        for j_depth in range(i_depth,i_depth+2):                                            # 右侧,上底
                            x.append(plot_depth[i_point+1,j_depth]-plot_depth[i_point+1,i_depth])
                            y.append(plot_grid[i_point+1]-plot_grid[i_point])
                            z.append(normal_velocities[i_point+1,j_depth])
                        for j_depth in range(i_depth,i_depth_temp+1):                                       # 左侧,下底
                            x.append(plot_depth[i_point,j_depth]-plot_depth[i_point,i_depth])
                            y.append(0)
                            z.append(normal_velocities[i_point,j_depth])
                        x0 = (2*a+b)*h/3/(a+b)
                        y0 = (a**2+b**2+a*b)/3/(a+b)
                        # 插值到形心位置
                        valid_speed[i_point, i_depth] = interp_IDW(x,y,z,x0,y0)
                        valid_areas[i_point, i_depth] = ((bottom_left - top_left) + (bottom_right - top_right)) * point_widths_m[i_point] / 2  # m²
                    break
                else:  # 不是底部,速度乘矩形面积
                    valid_areas[i_point, i_depth] = point_widths_m[i_point]*(bottom_left-top_left)          # m²
                    valid_speed[i_point, i_depth] = (normal_velocities[i_point, i_depth]+normal_velocities[i_point+1, i_depth]+normal_velocities[i_point+1, i_depth]+normal_velocities[i_point+1, i_depth+1])/4
                i_depth +=1

        # 验证有效面积
        max_area = np.nanmax(valid_areas)
        min_area = np.nanmin(valid_areas[valid_areas > 0])
        print(f"最大有效面积: {max_area:.2f} m²")
        print(f"最小有效面积: {min_area:.2f} m²")

        # === 7. 计算海水通量 ===
        print("\n计算海水通量...")

        # 计算通量 (Sv)
        total_flux_sv = np.nansum(valid_speed * valid_areas)/ 1e6

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
        plt.fill_between(plot_grid/1000., seafloor_depth, max_y,
                         color='lightgray', zorder=10)

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
        plt.savefig('output/fvcom/'+current_time_str+'normal_velocity_section.png', dpi=300)
        plt.show()

        # 存储通量数据入数组
        time_all_combined.append(current_time_str)
        flux_all_combined.append(total_flux_sv)

    with open('flux.txt', 'w') as fout:
        for i in range(len(time_all_combined)):
            fout.write(time_all_combined[i]+' '+str(flux_all_combined[i]))
            fout.write('\n')
    fout.close()