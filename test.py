import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dataframe = pd.read_csv(r"C:\Users\wsyzb\OneDrive\Desktop\CR1000XSeries_MIN_1.csv", header=0, usecols=[0,5,6], dtype=str)
dataframe['latitude'] = np.nan
dataframe['longitude'] = np.nan
for i_row in range(dataframe.shape[0]):
    dataframe.loc[i_row, 'Time'] = dataframe.loc[i_row, 'Time'][:-3]
    dataframe.loc[i_row, 'wind speed'] = float(dataframe.loc[i_row, 'wind speed'])
    dataframe.loc[i_row, 'wind direction'] = float(dataframe.loc[i_row, 'wind direction'])
print('read_csv ok')

with open(r"C:\Users\wsyzb\OneDrive\Desktop\tracks.gpx") as fin:
    lines = fin.readlines()
    i_line = 0
    while i_line<len(lines):
        if 'lat' in lines[i_line] and 'lon' in lines[i_line]:
            latitude = lines[i_line].split('"')[1]
            longitude = lines[i_line].split('"')[3]
            # 寻找时间
            i_line+=1
            print(i_line)
            time_str = lines[i_line][14:24]+' '+lines[i_line][25:30]
            for i_row in range(dataframe.shape[0]):
                if dataframe.loc[i_row, 'Time'] == time_str:
                    # 填入经纬度信息
                    dataframe.loc[i_row, 'latitude'] = float(latitude)
                    dataframe.loc[i_row, 'longitude'] = float(longitude)
                    break
            i_line+=1
        else:
            i_line+=1

dataframe['u'] = np.nan
dataframe['v'] = np.nan
for i_row in range(dataframe.shape[0]):
    dataframe.loc[i_row, 'u'] = float(dataframe.loc[i_row, 'wind speed'])*float(np.sin(np.deg2rad(dataframe.loc[i_row, 'wind direction'])))
    dataframe.loc[i_row, 'v'] = float(dataframe.loc[i_row, 'wind speed'])*float(np.cos(np.deg2rad(dataframe.loc[i_row, 'wind direction'])))

plt.figure(figsize=(10,10),dpi=100)
plt.quiver(dataframe['longitude'][::60],dataframe['latitude'][::60],dataframe['u'][::60],dataframe['v'][::60], scale=50)
plt.show()