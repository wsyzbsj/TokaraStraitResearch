import toml
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
plt.rcParams['font.family'] = 'Times New Roman, SimHei'
plt.rcParams['mathtext.fontset'] = 'stix'

# 主程序
if __name__ == "__main__":
    config = toml.load('config/config.toml')
    lat_min_zoom = config['zone_config']['lat_min_zoom']
    lat_max_zoom = config['zone_config']['lat_max_zoom']
    lon_min_zoom = config['zone_config']['lon_min_zoom']
    lon_max_zoom = config['zone_config']['lon_max_zoom']

    dataset = np.load('output/HYCOM/三维数据/all.npz')
    lat = dataset['lat']
    lon = dataset['lon']
    time = dataset['time']
    depth = dataset['depth']
    t = dataset['t']
    s = dataset['s']
    u = dataset['u']
    v = dataset['v']

    t_spring = []
    t_summer = []
    t_autumn = []
    t_winter = []

    s_spring = []
    s_summer = []
    s_autumn = []
    s_winter = []

    u_spring = []
    u_summer = []
    u_autumn = []
    u_winter = []

    v_spring = []
    v_summer = []
    v_autumn = []
    v_winter = []

    # 处理缺测
    t[t < -10] = np.nan
    s[s < -10] = np.nan
    u[u < -10] = np.nan
    v[v < -10] = np.nan

    # 按季节平均
    i_spring = 0
    i_summer = 0
    i_autumn = 0
    i_winter = 0
    for i_time,current_time in enumerate(time):
        if '12月' in current_time or '01月' in current_time or '02月' in current_time:     # 冬天
            t_winter.append(t[i_time])
            s_winter.append(s[i_time])
            u_winter.append(u[i_time])
            v_winter.append(v[i_time])
        elif '03月' in current_time or '04月' in current_time or '05月' in current_time:   # 春天
            t_spring.append(t[i_time])
            s_spring.append(s[i_time])
            u_spring.append(u[i_time])
            v_spring.append(v[i_time])
        elif '06月' in current_time or '07月' in current_time or '08月' in current_time:   # 夏天
            t_summer.append(t[i_time])
            s_summer.append(s[i_time])
            u_summer.append(u[i_time])
            v_summer.append(v[i_time])
        elif '09月' in current_time or '10月' in current_time or '11月' in current_time:   # 秋天
            t_autumn.append(t[i_time])
            s_autumn.append(s[i_time])
            u_autumn.append(u[i_time])
            v_autumn.append(v[i_time])

    t_winter = np.array(t_winter)
    s_winter = np.array(s_winter)
    u_winter = np.array(u_winter)
    v_winter = np.array(v_winter)
    t_spring = np.array(t_spring)
    s_spring = np.array(s_spring)
    u_spring = np.array(u_spring)
    v_spring = np.array(v_spring)
    t_summer = np.array(t_summer)
    s_summer = np.array(s_summer)
    u_summer = np.array(u_summer)
    v_summer = np.array(v_summer)
    t_autumn = np.array(t_autumn)
    s_autumn = np.array(s_autumn)
    u_autumn = np.array(u_autumn)
    v_autumn = np.array(v_autumn)

    ave_u_spring = np.nanmean(u_spring[:],axis=0)
    ave_u_summer = np.nanmean(u_summer[:],axis=0)
    ave_u_autumn = np.nanmean(u_autumn[:],axis=0)
    ave_u_winter = np.nanmean(u_winter[:],axis=0)

    ave_v_spring = np.nanmean(v_spring[:],axis=0)
    ave_v_summer = np.nanmean(v_summer[:],axis=0)
    ave_v_autumn = np.nanmean(v_autumn[:],axis=0)
    ave_v_winter = np.nanmean(v_winter[:],axis=0)

    ave_t_spring = np.nanmean(t_spring[:],axis=0)
    ave_t_summer = np.nanmean(t_summer[:],axis=0)
    ave_t_autumn = np.nanmean(t_autumn[:],axis=0)
    ave_t_winter = np.nanmean(t_winter[:],axis=0)

    ave_s_spring = np.nanmean(s_spring[:],axis=0)
    ave_s_summer = np.nanmean(s_summer[:],axis=0)
    ave_s_autumn = np.nanmean(s_autumn[:],axis=0)
    ave_s_winter = np.nanmean(s_winter[:],axis=0)

    # 绘制SST场与速度场
    for i_depth in range(len(depth)):
        # 创建图形并调整布局，为底部colorbar腾出空间
        fig = plt.figure(figsize=(10, 10))
        plt.subplots_adjust(bottom=0.08, top=0.90, right=0.85)  # 调整底部边距
        plt.suptitle('Water Temperature and Speed Field at '+str(int(depth[i_depth]))+' m depth', fontsize=20)

        # 获取水温最大最小值
        vmax = np.ceil(np.nanmax([np.nanmax(ave_t_spring[i_depth,:,:]),np.nanmax(ave_t_summer[i_depth,:,:]),np.nanmax(ave_t_autumn[i_depth,:,:]),np.nanmax(ave_t_winter[i_depth,:,:])]))
        vmin = np.floor(np.nanmin([np.nanmin(ave_t_spring[i_depth,:,:]),np.nanmin(ave_t_summer[i_depth,:,:]),np.nanmin(ave_t_autumn[i_depth,:,:]),np.nanmin(ave_t_winter[i_depth,:,:])]))

        # 春季
        ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
        ax1.set_extent([lon_min_zoom,lon_max_zoom,lat_min_zoom,lat_max_zoom], crs=ccrs.PlateCarree())
        # 创建正确的二维网格坐标
        X_grid, Y_grid = np.meshgrid(lon, lat, indexing='ij')

        # 创建颜色映射
        cmap = cmaps.MPL_RdBu[::-1]
        cmap.set_bad('gray', 1.0)

        # 绘制水温数据 (使用二维坐标)
        mesh = ax1.contourf(
            X_grid, Y_grid, ave_t_spring[i_depth,:,:],  # 使用二维坐标和流速数据
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=np.arange(vmin, vmax, 0.05),
            zorder=0
        )

        # 添加等值线
        contour = ax1.contour(
            X_grid, Y_grid, ave_t_spring[i_depth,:,:],
            levels=np.arange(-6,40,2),
            colors='red',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 添加陆地地形
        ax1.add_feature(cfeature.LAND, color='darkgray')
        ax1.clabel(contour,inline=True,fontsize=10)
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = True

        # 设置坐标轴
        ax1.set_ylabel(r'Latitude $(\degree)$', fontsize=12)
        ax1.set_title('Spring', fontsize=12)

        # 夏季
        ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
        ax2.set_extent([lon_min_zoom,lon_max_zoom,lat_min_zoom,lat_max_zoom], crs=ccrs.PlateCarree())

        # 绘制水温数据 (使用二维坐标)
        mesh = ax2.contourf(
            X_grid, Y_grid, ave_t_summer[i_depth,:,:],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=np.arange(vmin, vmax, 0.05),
            zorder=0
        )

        # 添加等值线
        contour = ax2.contour(
            X_grid, Y_grid, ave_t_summer[i_depth,:,:],
            levels=np.arange(-6,40,2),
            colors='red',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 添加陆地
        ax2.add_feature(cfeature.LAND, color='darkgray')
        ax2.clabel(contour,inline=True,fontsize=10)
        gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = False

        # 设置坐标轴
        ax2.set_title('Summer', fontsize=12)

        # 秋季
        ax3 = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
        ax3.set_extent([lon_min_zoom,lon_max_zoom,lat_min_zoom,lat_max_zoom], crs=ccrs.PlateCarree())

        # 绘制水温数据 (使用二维坐标)
        mesh = ax3.contourf(
            X_grid, Y_grid, ave_t_autumn[i_depth,:,:],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=np.arange(vmin, vmax, 0.05),
            zorder=0
        )

        # 添加等值线
        contour = ax3.contour(
            X_grid, Y_grid, ave_t_autumn[i_depth,:,:],
            levels=np.arange(-6,40,2),
            colors='red',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 添加陆地
        ax3.add_feature(cfeature.LAND, color='darkgray')
        ax3.clabel(contour,inline=True,fontsize=10)
        gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = True

        # 设置坐标轴
        ax3.set_xlabel(r'Longitude $(\degree)$', fontsize=12)
        ax3.set_ylabel(r'Latitude $(\degree)$', fontsize=12)
        ax3.set_title('Autumn', fontsize=12)

        # 冬季
        ax4 = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
        ax4.set_extent([lon_min_zoom,lon_max_zoom,lat_min_zoom,lat_max_zoom], crs=ccrs.PlateCarree())

        # 绘制水温数据 (使用二维坐标)
        mesh = ax4.contourf(
            X_grid, Y_grid, ave_t_winter[i_depth,:,:],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            levels=np.arange(vmin, vmax, 0.05),
            zorder=0
        )

        # 添加等值线
        contour = ax4.contour(
            X_grid, Y_grid, ave_t_winter[i_depth,:,:],
            levels=np.arange(-6,40,2),
            colors='red',
            linewidths=0.5,
            linestyles='--',
            zorder=0
        )

        # 添加陆地
        ax4.add_feature(cfeature.LAND, color='darkgray')
        ax4.clabel(contour,inline=True,fontsize=10)
        gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = False

        # 设置坐标轴
        ax4.set_xlabel(r'Longitude $(\degree)$', fontsize=12)
        ax4.set_title('Winter', fontsize=12)

        # 在右侧添加一个共享的colorbar（垂直方向）
        cbar_ax = fig.add_axes([0.88, 0.10, 0.02, 0.78])  # [left, bottom, width, height]
        cbar = fig.colorbar(mesh, cax=cbar_ax, ticks=np.arange(-6,40,2), orientation='vertical')
        cbar.set_label(r'Water temperature $(\degree C)$', fontsize=12)

        plt.savefig('output/HYCOM/平面/'+str(int(depth[i_depth]))+'m深度季节SST与流速.png', dpi=600)
        plt.show()
        plt.close()