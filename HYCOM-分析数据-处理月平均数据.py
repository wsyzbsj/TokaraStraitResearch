import numpy as np
from hycom.info import read_field_names
from hycom.io import read_hycom_fields
import netCDF4
import matplotlib.pyplot as plt

fields = read_field_names('/home/yzbsj/Data/海洋数据/HYCOM/expt_93.0_meanstd/regional.grid.a')
print(fields)
coordinates = read_hycom_fields('/home/yzbsj/Data/海洋数据/HYCOM/expt_93.0_meanstd/regional.grid.a', fields)
print(coordinates)
# fields = read_field_names('/home/yzbsj/Data/海洋数据/HYCOM/expt_93.0_meanstd/930_archMN.2023_01.b')
# dataset = netCDF4.Dataset('/home/yzbsj/Data/海洋数据/HYCOM/expt_93.0_meanstd/depth_GLBb0.08_09m11.nc', 'r')
# dep = dataset.variables['depth'][0][:]
# lat = dataset.variables['Latitude'][:]
# lon = dataset.variables['Longitude'][:]
# # 需要流速数据,u-vel和v-vel
# data = read_hycom_fields('/home/yzbsj/Data/海洋数据/HYCOM/expt_93.0_meanstd/930_archMN.2023_01.a', fields)
# u_vel = data['u-vel.']  # (41,3298,4500) 41层×3298纬度×4500经度
# v_vel = data['v-vel.']  # (41,3298,4500) 41层×3298纬度×4500经度
# plt.contourf(lon,lat,(u_vel[0]**2+v_vel[0]**2)**0.5,cmap=plt.get_cmap('jet'))
# plt.colorbar()
# plt.show()
print(1)