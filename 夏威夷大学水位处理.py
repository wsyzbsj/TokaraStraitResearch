import numpy as np
import pandas as pd
import glob
import pathlib
import netCDF4
from datetime import datetime,timedelta
from functools import reduce
from geopy.distance import geodesic as gd
from scipy.interpolate import RegularGridInterpolator

def calc_distance(lat1, lon1, lat2, lon2):  # 计算两点之间的距离
    # 两点的经度和纬度
    p1 = (lat1, lon1)
    p2 = (lat2, lon2)
    d = gd(p1, p2).m	# 得到米为单位的结果
    return d

station_used = ['Naze', 'Nakano Shima','Nishinoomote']
for i_station in range(len(station_used)):
    exec("dataframe_"+str(i_station)+" = pd.DataFrame()")
    # exec("dataframe_"+str(i_station)+"['station'] = ''")
    # exec("dataframe_"+str(i_station)+"['time'] = ''")
    # exec("dataframe_"+str(i_station)+"['latitude'] = ''")
    # exec("dataframe_"+str(i_station)+"['longitude'] = ''")
    # exec("dataframe_"+str(i_station)+"['sea level'] = ''")

files = sorted(glob.glob("/home/yzbsj/Data/海洋数据/夏威夷大学验潮站数据/*.nc"))
lat_sta = np.zeros((len(station_used)),dtype=float)    # 站点纬度
lon_sta = np.zeros((len(station_used)),dtype=float)    # 站点经度
for file in files:
    station = pathlib.Path(file).stem
    i_sta = np.nan
    if station in station_used:
        for i_station in range(len(station_used)):
            if station_used[i_station] == station:
                i_sta = i_station
                break
        dataset = netCDF4.Dataset(file, 'r')
        time = dataset.variables['time'][:]
        lat = dataset.variables['lat'][:]
        if len(lat) == 1:
            lat_sta[i_sta] = lat[0]
        lon = dataset.variables['lon'][:]
        if len(lon) == 1:
            lon_sta[i_sta] = lon[0]
        sea_level = dataset.variables['sea_level'][0,:]
        time_start = datetime.strptime(dataset.variables['time'].units[11:], "%Y-%m-%d %H:%M:%S")
        real_time = []
        for i_time in range(len(time)):
            real_time.append(time_start + timedelta(days=time[i_time]))
        exec("dataframe_"+str(i_sta)+".loc[:,'time'] = real_time")
        exec("dataframe_"+str(i_sta)+".loc[:,'"+station+" sea level'] = sea_level")
        exec("dataframe_"+str(i_sta)+".set_index('time')")
    else:
        continue

dataframes = [dataframe_0, dataframe_1, dataframe_2]
merged_df = reduce(lambda left, right: pd.merge(left, right, on='time', how='outer'), dataframes)
merged_df.to_excel("output/验潮站数据/合并.xlsx")

# 计算地转输运
dataset = netCDF4.Dataset('/home/yzbsj/Data/其他数据/GEBCO/gebco_2024/GEBCO_2024.nc', 'r')
depth_lon = dataset.variables['lon'][:]
depth_lat = dataset.variables['lat'][:]
depth_elevation = -dataset.variables['elevation'][:]
dataset.close()
num_points = 100
for i in range(len(station_used)-1):
    distance = calc_distance(lat_sta[i], lon_sta[i], lat_sta[i+1], lon_sta[i+1])
    lat_list = np.linspace(lat_sta[i], lat_sta[i+1], num=num_points)
    lon_list = np.linspace(lon_sta[i], lon_sta[i+1], num=num_points)
    points = np.column_stack((lat_list, lon_list))
    # 创建水深插值函数并插值
    depth_interp = RegularGridInterpolator(
        (depth_lat, depth_lon),
        depth_elevation,
        method='linear',
        bounds_error=False,
        fill_value=np.nan
    )
    seafloor_depth = depth_interp(points)
    average_depth = np.nanmean(seafloor_depth)  # 两站之间的平均深度

    time = merged_df.shape[0]  # 总时间长度

    # 生成模拟潮位数据 (两个站有相位差)
    eta_A = merged_df['Naze sea level']
    eta_B = merged_df['Nakano Shima sea level']

    # KVT10=0.36×SLDA10(石垣-基隆)+20.89
    flow = 0.36*(eta_B-eta_A)+20.89
    print(1)