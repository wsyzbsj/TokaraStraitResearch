import pandas as pd
import glob
import pathlib
import netCDF4
from datetime import datetime,timedelta

station_used = ['Ishigaki','Keelung']
for i_station in range(len(station_used)):
    exec("dataframe_"+str(i_station)+" = pd.DataFrame()")
    # exec("dataframe_"+str(i_station)+"['station'] = ''")
    # exec("dataframe_"+str(i_station)+"['time'] = ''")
    # exec("dataframe_"+str(i_station)+"['latitude'] = ''")
    # exec("dataframe_"+str(i_station)+"['longitude'] = ''")
    # exec("dataframe_"+str(i_station)+"['sea level'] = ''")

files = sorted(glob.glob("/home/yzbsj/Data/海洋数据/夏威夷大学验潮站数据/*.nc"))
for file in files:
    station = pathlib.Path(file).stem
    if station in station_used:
        for i_station in range(len(station_used)):
            if station_used[i_station] == station:
                i_sta = i_station
                break
        dataset = netCDF4.Dataset(file, 'r')
        time = dataset.variables['time'][:]
        lat = dataset.variables['lat'][:]
        if len(lat) == 1:
            lat = lat[0]
        lon = dataset.variables['lon'][:]
        if len(lon) == 1:
            lon = lon[0]
        sea_level = dataset.variables['sea_level'][0,:]
        time_start = datetime.strptime(dataset.variables['time'].units[11:], "%Y-%m-%d %H:%M:%S")
        real_time = []
        for i_time in range(len(time)):
            real_time.append(time_start + timedelta(days=time[i_time]))
        exec("dataframe_"+str(i_sta)+".loc[:,'time'] = real_time")
        exec("dataframe_"+str(i_sta)+".loc[:,'"+station+"sea level'] = sea_level")
        exec("dataframe_"+str(i_sta)+".loc[:,'"+station+"station'] = station")
        exec("dataframe_"+str(i_sta)+".loc[:,'"+station+"latitude'] = lat")
        exec("dataframe_"+str(i_sta)+".loc[:,'"+station+"longitude'] = lon")
        exec("dataframe_"+str(i_sta)+".set_index('time')")
    else:
        continue

combined_df = pd.merge(dataframe_0, dataframe_1, on='time', how='outer')
combined_df.to_excel("output/验潮站数据/合并.xlsx")