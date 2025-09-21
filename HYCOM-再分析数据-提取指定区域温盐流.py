import glob
import pathlib
import os
import netCDF4
from datetime import datetime,timedelta
import numpy as np
import toml

# 主程序
if __name__ == "__main__":
    speed_all_combined = []
    flux_all_combined = []
    time_all_combined = []

    # 读取设置——注:用不到大区域故将lat_max等变量名给数据范围
    config = toml.load('config/config.toml')

    lat_min_zoom = config['zone_config']['lat_min_zoom']
    lat_max_zoom = config['zone_config']['lat_max_zoom']
    lon_min_zoom = config['zone_config']['lon_min_zoom']
    lon_max_zoom = config['zone_config']['lon_max_zoom']
    lat_max = lat_max_zoom+0.5
    lat_min = lat_min_zoom-0.5
    lon_max = lon_max_zoom+0.5
    lon_min = lon_min_zoom-0.5

    files_uv = sorted(glob.glob(config['file_config']['HYCOM_reanalysis_uv']))
    files_ts = sorted(glob.glob(config['file_config']['HYCOM_reanalysis_ts']))

    u_all = []
    v_all = []
    t_all = []
    s_all = []
    lat_all = []
    lon_all = []
    time_all = []
    depth_all = []

    for i_file,file in enumerate(files_uv):
        filename_uv = pathlib.Path(file).name
        file_uv = os.path.join(pathlib.Path(file).parent, filename_uv)
        filename_ts = filename_uv.replace('uv','ts')
        file_ts = os.path.join(pathlib.Path(file).parent, filename_ts)
        dataset_uv = netCDF4.Dataset(file_uv,'r')
        dataset_ts = netCDF4.Dataset(file_ts,'r')

        uv_lat = dataset_uv.variables['lat'][:]
        uv_lon = dataset_uv.variables['lon'][:]
        lat_idx_min = np.argwhere(uv_lat > lat_min-1)[0][0]
        lat_idx_max = np.argwhere(uv_lat < lat_max+1)[-1][0]
        lon_idx_min = np.argwhere(uv_lon > lon_min-1)[0][0]
        lon_idx_max = np.argwhere(uv_lon < lon_max+1)[-1][0]
        uv_lat = uv_lat[lat_idx_min:lat_idx_max]
        uv_lon = uv_lon[lon_idx_min:lon_idx_max]
        uv_depth = dataset_uv.variables['depth'][:]
        uv_time = pathlib.Path(file_uv).name.split('.')[2].split('_')[0]+'年'+pathlib.Path(file_uv).name.split('.')[2].split('_')[1]+'月'
        u = dataset_uv.variables['water_u'][0, :, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]
        v = dataset_uv.variables['water_v'][0, :, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]

        ts_lat = dataset_ts.variables['lat'][:]
        ts_lon = dataset_ts.variables['lon'][:]
        lat_idx_min = np.argwhere(ts_lat > lat_min-1)[0][0]
        lat_idx_max = np.argwhere(ts_lat < lat_max+1)[-1][0]
        lon_idx_min = np.argwhere(ts_lon > lon_min-1)[0][0]
        lon_idx_max = np.argwhere(ts_lon < lon_max+1)[-1][0]
        ts_lat = ts_lat[lat_idx_min:lat_idx_max]
        ts_lon = ts_lon[lon_idx_min:lon_idx_max]
        ts_depth = dataset_ts.variables['depth'][:]
        ts_time = pathlib.Path(file_ts).name.split('.')[2].split('_')[0]+'年'+pathlib.Path(file_ts).name.split('.')[2].split('_')[1]+'月'
        t = dataset_ts.variables['water_temp'][0, :, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]
        s = dataset_ts.variables['salinity'][0, :, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max]

        if i_file == 0:
            depth_all = ts_depth
            for i in range(len(uv_lat)):
                lat_all.append(round(uv_lat[i],2))
                lon_all.append(round(uv_lon[i],2))
        u_all.append(u)
        v_all.append(v)
        t_all.append(t)
        s_all.append(s)
        time_all.append(uv_time)

    lat_all = np.array(lat_all)
    lon_all = np.array(lon_all)
    time_all = np.array(time_all)
    depth_all = np.array(depth_all)
    u_all = np.array(u_all)
    v_all = np.array(v_all)
    t_all = np.array(t_all)
    s_all = np.array(s_all)
    np.savez('output/HYCOM/三维数据/all',lat=lat_all,lon=lon_all,depth=depth_all,t=t_all,s=s_all,u=u_all,v=v_all,time=time_all)