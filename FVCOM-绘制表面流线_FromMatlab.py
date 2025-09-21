import toml
import netCDF4
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy.interpolate import griddata
from datetime import datetime, timedelta
import pandas as pd
import cartopy.io.shapereader as shpreader
from shapely.geometry import box
from shapely.prepared import prep
from shapely.strtree import STRtree
import matplotlib.path as mpath
from matplotlib.patches import Rectangle
from scipy.integrate import odeint
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import os

# 确保输出目录存在
os.makedirs("output/FVCOM", exist_ok=True)

# 设置
config = toml.load('config/config.toml')
plt.rcParams['font.family'] = 'Times New Roman, Microsoft YaHei'
plt.rcParams['mathtext.fontset'] = 'stix'
# 局部放大区域
lon_min = config['zone_config']['lon_min_zoom']
lon_max = config['zone_config']['lon_max_zoom']
lat_min = config['zone_config']['lat_min_zoom']
lat_max = config['zone_config']['lat_max_zoom']
proj = ccrs.PlateCarree()

# 将经纬度转换为投影坐标（等效于 m_ll2xy）
def lonlat_to_xy(lon, lat, projection=proj):
    """将经纬度坐标转换为投影坐标"""
    return projection.transform_points(ccrs.PlateCarree(), np.array(lon), np.array(lat))

def id2axis(distance,id):
    n=int(np.ceil((0.5-distance/2)/distance)*2+1)        # 分割的数量
    min_distance=(1-(n-2)*distance)/2               # 两端最小的距离
    if id==1:
        xpoint=min_distance/2
    elif id==n:
        xpoint=1-min_distance/2
    else:
        xpoint=min_distance+(id-1.5)*distance
    return xpoint

def sub2ind(sizeArray, *subs):
    """
    模拟 MATLAB 的 sub2ind 函数

    参数:
        sizeArray: 数组的形状
        *subs: 下标索引（可以是多个参数，每个参数是一个下标数组）

    返回:
        线性索引（MATLAB 风格，从1开始）
    """
    # 将下标转换为 0-based
    subs_0based = [s - 1 for s in subs]

    # 计算线性索引
    ind = np.ravel_multi_index(subs_0based, dims=sizeArray, order='F')

    # 转换为 1-based（MATLAB 风格）
    return ind + 1

def axis2id(x,y,distance):
    N=int(np.ceil((0.5-distance/2)/distance)*2+1)                                                                       # 分割的数量
    min_distance=(1-(N-2)*distance)/2                                                                                   # 两端最小的距离

    # x的位置
    if x<=min_distance:
        pos_id_x=1
    elif x>=1-min_distance:
        pos_id_x=N
    else:
        pos_id_x=int(np.ceil((x-min_distance)/distance)+1)

    # y的位置
    if y<=min_distance:
        pos_id_y=1
    elif y>=1-min_distance:
        pos_id_y=N
    else:
        pos_id_y=int(np.ceil((y-min_distance)/distance)+1)

    # xy转ind
    pos_id=sub2ind([N,N],pos_id_y,pos_id_x)

    return pos_id

def delete_self(sl_i, xy_end, dend, xy_start, dstart):
    """
    删除自相交或过于接近的流线点，并更新网格状态

    参数:
        sl_i: 流线数据 (N行2列的数组，包含x,y坐标)
        xy_end: 终点网格
        dend: 终点网格间距
        xy_start: 起点网格
        dstart: 起点网格间距

    返回:
        sl_i: 处理后的流线
        xy_end: 更新后的终点网格
        xy_start: 更新后的起点网格
    """
    # 删除包含 NaN 的行
    del_in = np.isnan(sl_i[:, 0])
    sl_i = sl_i[~del_in, :]

    # 获取流线点数
    N = sl_i.shape[0]
    if N == 0:
        print('流线点数为0')
        exit(-1)
        return sl_i, xy_end, xy_start

    # 获取第一个点的网格索引并标记
    pos_id_last = axis2id(sl_i[0, 0], sl_i[0, 1], dend)
    xy_end[pos_id_last] = 1  # 标记第一个点

    # 顺便标记起点网格
    pos_id_s = axis2id(sl_i[0, 0], sl_i[0, 1], dstart)
    xy_start[pos_id_s] = 1

    # 遍历流线上的其他点
    j = 1
    while j < N:
        pos_id_now = axis2id(sl_i[j, 0], sl_i[j, 1], dend)
        if pos_id_now != pos_id_last:
            # 如果现在的点和原有的点在同一区域，则不管它
            # 如果不在同一区域，检测新的点是否已经被占用
            if xy_end[pos_id_now] == 1:     # 如果该点被占用，说明出现与其它流线太近的情况，则直接停止
                j -= 1
                break
            else:                           # 如果没被占用，则把新点添加上
                xy_end[pos_id_now] = 1
                pos_id_last = pos_id_now
        # 顺便标记起点网格
        pos_id_s = axis2id(sl_i[j, 0], sl_i[j, 1], dstart)
        xy_start[pos_id_s] = 1
        j += 1
    # 删除j之后的所有点
    if j < N:
        sl_i = sl_i[:j+1, :]
    return sl_i, xy_end, xy_start

def matlab_histcounts_indices(data, bins):
    # 使用 digitize 获取 bin 索引
    indices = np.digitize(data, bins)

    # 处理超出范围的值（MATLAB 中返回 0）
    # 对于小于最小 bin 的值，digitize 返回 0，与 MATLAB 一致
    # 对于大于最大 bin 的值，digitize 返回 len(bins)+1，需要调整为 0
    indices[indices > len(bins)] = 0

    return indices

def my_streamline_mutli(x,y,u,v,dstart:float,num:int,magnify:int) -> list:
    streamline_sum = []
    dend = 0.5*dstart
    xmin=x.min()
    xmax=x.max()
    ymin=y.min()
    ymax=y.max()

    # 归一化，将流场缩放为0-1区间的矩形
    xn=(x-xmin)/(xmax-xmin)
    yn=(y-ymin)/(ymax-ymin)
    un=u/(xmax-xmin)
    vn=v/(ymax-ymin)

    num_start=int(np.ceil((0.5-dstart/2)/dstart)*2+1)
    num_end=int(np.ceil((0.5-dend/2)/dend)*2+1)

    # 初始化所有网格点，0代表可以放置新点，1代表已经存在原有的点
    xy_start = np.zeros((num_start,num_start),dtype=float)
    xy_end = np.zeros((num_end,num_end),dtype=float)

    # 标记陆地上点为1:，代表不会进入循环
    mask = np.isnan(un).astype(float)
    zy = np.linspace(0,1,xy_start.shape[0])
    zx = np.linspace(0,1,xy_start.shape[1])
    zxx,zyy = np.meshgrid(zx,zy)
    # 原本插值:xy_start = interp2(xn,yn,mask,zxx,zyy)
    # 创建插值器
    interpolator = RegularGridInterpolator((yn[0, :], xn[:, 0]), mask, method='linear', bounds_error=False, fill_value=np.nan)
    # 进行插值
    xy_start = interpolator((zyy, zxx))  # 注意坐标顺序

    py_in = np.where(xy_start.flatten(order='C')>0)
    for i_py_in in py_in:
        i_x = i_py_in//xy_start.shape[1]
        i_y = i_py_in%xy_start.shape[1]
        xy_start[i_x,i_y] = 1
    zy = np.linspace(0,1,xy_end.shape[0])
    zx = np.linspace(0,1,xy_end.shape[1])
    zxx,zyy = np.meshgrid(zx,zy)
    # 原本的插值:xy_end = interp2(xn,yn,mask,zxx,zyy)
    # 创建插值器
    interpolator = RegularGridInterpolator((yn[0, :], xn[:, 0]), mask, method='linear', bounds_error=False, fill_value=np.nan)
    # 进行插值
    xy_end = interpolator((zyy, zxx))  # 注意坐标顺序
    py_in = np.where(xy_end.flatten(order='C')>0)
    for i_py_in in py_in:
        i_x = i_py_in//xy_end.shape[1]
        i_y = i_py_in%xy_end.shape[1]
        xy_end[i_x,i_y] = 1

    # 将流线划分为num种，速度越大的流线越长
    length_sl=np.linspace(5,40,num)                                                                                     # 按速度分类
    V2=(un**2+vn**2)**0.5
    V2_max=V2.max()
    V2_min=V2.min()

    V2_space=np.linspace(V2_min,V2_max,num+1)

    # 1. 当xy_start内还有可放置的新点的位置时，循环
    k=0                                                                                                                 # 流线数(循环次数)
    streamline_seed = []
    while not np.all(xy_start):
        # 2. 随机一个start内网格点作为种子点
        [start_id_y,start_id_x]=np.where(xy_start==0)
        randnum=np.random.randint(1, start_id_y.shape[0]+1)
        x_pos_i=id2axis(dstart,start_id_x[randnum])
        y_pos_i=id2axis(dstart,start_id_y[randnum])
        streamline_seed.append([x_pos_i,y_pos_i])                                                                       # 保存种子点
        # 原本的插值:V2_seed=interp2(xn,yn,V2,x_pos_i,y_pos_i)                                                               # 计算种子点处的速度
        # 创建插值器
        interpolator = RegularGridInterpolator((yn[0, :], xn[:, 0]), V2, method='linear', bounds_error=False, fill_value=np.nan)
        # 进行插值
        V2_seed = interpolator((zyy, zxx))  # 注意坐标顺序
        sl_N = matlab_histcounts_indices(V2_seed, V2_space)
        # [~,~,sl_N] = histcounts(V2_seed,V2_space)
        for i_x in range(sl_N.shape[0]):
            for i_y in range(sl_N.shape[1]):
                if sl_N[i_x,i_y]==0 or np.isnan(sl_N[i_x,i_y]):
                    sl_N[i_x,i_y] = 5

        num_streamline=round(length_sl[sl_N]) * magnify                                                                 # 流线总数量

        # 3. 绘制流线
        # MATLAB原代码
        # streamline_i_1 = stream2(xn,yn, un, vn,x_pos_i,y_pos_i,[0.1,num_streamline])
        # streamline_i_2 = stream2(xn,yn,-un,-vn,x_pos_i,y_pos_i,[0.1,num_streamline])
        def velocity_func(r, t, interpolator_u, interpolator_v):                                                        # 定义速度场函数
            x, y = r
            # 使用插值获取速度
            u_val = interpolator_u([x, y])[0]
            v_val = interpolator_v([x, y])[0]
            return [u_val, v_val]

        # 创建速度场插值器
        interpolator_u = RegularGridInterpolator((yn[0,:], xn[:,0]), un.T)
        interpolator_v = RegularGridInterpolator((yn[0,:], xn[:,0]), vn.T)

        # 计算正向流线（类似 streamline_i_1）
        t1 = np.linspace(0, num_streamline/10, num_streamline*10)
        streamline_i_1 = odeint(velocity_func, [x_pos_i, y_pos_i], t1, args=(interpolator_u, interpolator_v))

        # 计算反向流线（类似 streamline_i_2）
        t2 = np.linspace(0, -num_streamline/10, num_streamline*10)
        streamline_i_2 = odeint(velocity_func, [x_pos_i, y_pos_i], t2, args=(interpolator_u, interpolator_v))
        # 4. 以xy_end为标准，删除自相交或间隔太近的点。并顺便标记xy_end
        [streamline_i_1,xy_end,xy_start]=delete_self(streamline_i_1[0], xy_end,dend,xy_start,dstart)
        [streamline_i_2,xy_end,xy_start]=delete_self(streamline_i_2[0], xy_end,dend,xy_start,dstart)
        # 5. 保存
        streamline_k=[streamline_i_2.flipud(),streamline_i_1[:,:]]                                                      # 新的流线
        streamline_sum.append([])
        streamline_sum[k]=[xmin+streamline_k[:][0]*(xmax-xmin), ymin+streamline_k[:][1]*(ymax-ymin)]                    # 从归一化还原
        k += 1
    streamline_seed = np.array(streamline_seed)
    streamline_sum = np.array(streamline_sum)
    streamline_seed=[streamline_seed[:,0]*(xmax-xmin)+xmin, streamline_seed[:,1]*(ymax-ymin)+ymin]
    return streamline_sum,streamline_seed

# 读取流场数据
dataset = netCDF4.Dataset(config['file_config']['fvcom_file'], 'r')
nv = dataset.variables['nv'][:]
time = dataset.variables['time'][0]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
latc = dataset.variables['latc'][:]
lonc = dataset.variables['lonc'][:]
u = dataset.variables['u'][:]
v = dataset.variables['v'][:]
times = dataset.variables['Times'][:]
h = dataset.variables['h'][:]  # 水深
zeta = dataset.variables['zeta'][:]  # 水位
siglay = dataset.variables['siglay'][:]  # sigma层
temp = dataset.variables['temp'][:,0,:]  # 表面温度

# 创建规则网格
print('创建网格')
grid_x, grid_y = np.mgrid[lon_min:lon_max:200j, lat_min:lat_max:200j]  # 分辨率

print('读取海岸线数据')
shp_path = config['file_config']['polygon_file']
reader = shpreader.Reader(shp_path)

# 创建空间索引
geoms = list(reader.geometries())
tree = STRtree(geoms)

# 创建目标区域边界框
target_bbox = box(lon_min, lat_min, lon_max, lat_max)
prep_bbox = prep(target_bbox)  # 预处理加速查询

# 查询相交的几何体
intersect_indices = tree.query(target_bbox, predicate="intersects")
valid_geoms = [geoms[i] for i in intersect_indices]

# 使用更高效的方法计算陆地掩膜
print('使用matplotlib.Path计算陆地掩膜')
land_mask = np.zeros(grid_x.shape, dtype=bool)  # 初始化陆地掩膜为False

# 收集所有多边形的顶点
all_polygons = []
for geom in valid_geoms:
    if geom.geom_type == 'Polygon':
        all_polygons.append(geom)
    elif geom.geom_type == 'MultiPolygon':
        all_polygons.extend(geom.geoms)

# 创建包含所有多边形的复合路径
vertices = []
codes = []
for poly in all_polygons:
    # 处理外部环
    ext_vertices = np.array(poly.exterior.coords)
    vertices.append(ext_vertices)
    codes.append([mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(ext_vertices) - 2) + [mpath.Path.CLOSEPOLY])

    # 处理内部环（孔洞）
    for interior in poly.interiors:
        int_vertices = np.array(interior.coords)
        vertices.append(int_vertices)
        codes.append([mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(int_vertices) - 2) + [mpath.Path.CLOSEPOLY])

# 创建复合路径
compound_path = mpath.Path(np.concatenate(vertices), np.concatenate(codes))

# 一次性检查所有网格点
points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
inside = compound_path.contains_points(points)
land_mask = inside.reshape(grid_x.shape)

# 分别获取目标区域内的节点索引和单元索引
print('获取目标区域内的节点和单元索引')
valid_node_indices = []  # 用于温度数据（节点数据）
valid_cell_indices = []  # 用于流速数据（单元数据）

# 获取有效节点索引（用于温度）
for i in range(len(lat)):
    if lon_min - 0.5 < lon[i] < lon_max + 0.5 and lat_min - 0.5 < lat[i] < lat_max + 0.5:
        valid_node_indices.append(i)
valid_node_indices = np.array(valid_node_indices)

# 获取有效单元索引（用于流速）
for i in range(len(latc)):
    if lon_min - 0.5 < lonc[i] < lon_max + 0.5 and lat_min - 0.5 < latc[i] < lat_max + 0.5:
        valid_cell_indices.append(i)
valid_cell_indices = np.array(valid_cell_indices)

# 处理时间数据
print('处理时间数据')
times_all = []
for i_time in range(len(times)):
    temp_str = ''.join(times[i_time, j].decode('utf-8') for j in range(len(times[0])))
    times_all.append((pd.to_datetime(temp_str) + timedelta(hours=8)).strftime('%Y-%m-%d %H时'))

# 设置投影和范围
print('设置投影与绘图范围')
extents = [lon_min, lon_max, lat_min, lat_max]

for i_time in range(u.shape[0]):
    print(f"正在处理时间步: {i_time + 1}/{u.shape[0]} - {times_all[i_time]}")

    # 创建图形 - 设置背景为灰色
    fig = plt.figure(figsize=(12, 10), dpi=150)

    ax = plt.axes(projection=proj)
    ax.set_extent(extents, crs=proj)
    ax.set_facecolor('lightgray')  # 设置坐标轴区域背景为灰色
    ax.set_title(times_all[i_time] + " 表面温度与流线", fontsize=16)

    # 获取当前时间步的数据
    current_u = u[i_time, 0, valid_cell_indices]  # 表面层u分量（单元数据）
    current_v = v[i_time, 0, valid_cell_indices]  # 表面层v分量（单元数据）
    current_temp = temp[i_time, valid_node_indices]  # 表面温度（节点数据）
    current_speed = np.sqrt(current_u ** 2 + current_v ** 2)  # 流速大小

    # 插值到规则网格
    # 温度插值（使用节点位置）
    points_temp = np.column_stack((lon[valid_node_indices], lat[valid_node_indices]))
    grid_temp = griddata(points_temp, current_temp, (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 流速插值（使用单元中心位置）
    points_flow = np.column_stack((lonc[valid_cell_indices], latc[valid_cell_indices]))
    grid_u = griddata(points_flow, current_u, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_v = griddata(points_flow, current_v, (grid_x, grid_y), method='linear', fill_value=np.nan)
    grid_speed = griddata(points_flow, current_speed, (grid_x, grid_y), method='linear', fill_value=np.nan)

    # 应用陆地掩膜
    grid_temp[land_mask] = np.nan
    grid_u[land_mask] = np.nan
    grid_v[land_mask] = np.nan
    grid_speed[land_mask] = np.nan

    # 绘制温度填色图 - 设置在最底层
    temp_cmap = plt.get_cmap('coolwarm')  # 使用适合温度的颜色映射
    cf_temp = ax.contourf(
        grid_x, grid_y, grid_temp,
        levels=50,
        cmap=temp_cmap,
        transform=proj,
        zorder=1  # 设置在最底层
    )

    # 添加温度colorbar
    cbar_temp = plt.colorbar(cf_temp, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
    cbar_temp.set_label('温度 (℃)', fontsize=12)

    # 只绘制相关几何体 - 移除填充，保留边界线
    for geom in valid_geoms:
        ax.add_geometries(
            [geom],
            crs=ccrs.PlateCarree(),
            facecolor='lightgray',  # 移除填充
            edgecolor='black',  # 保留黑色边界
            zorder=10
        )

    # 创建海洋掩膜（非陆地）
    ocean_mask = ~land_mask

    # 创建海洋区域的流速网格（陆地区域设为np.nan）
    grid_u_ocean = np.where(ocean_mask, grid_u, np.nan)
    grid_v_ocean = np.where(ocean_mask, grid_v, np.nan)

    # 确保grid_speed_smooth中没有NaN值
    grid_speed_filled = np.nan_to_num(grid_speed, nan=np.nan)

    streamline_sum,streamline_seed=my_streamline_mutli(grid_x,grid_y,grid_u_ocean,grid_v_ocean,0.04,10,5)

    xy_ratio = ax.get_data_ratio()
    xy_lim = ax.get_xlim(), ax.get_ylim()

    for j in range(streamline_sum.shape[1]):
        if streamline_sum[j].shape[0]<=1:
            continue        # 忽略异常流线
        ax.plot(streamline_sum[j][:,0],streamline_sum[j][:,1],'k','LineWidth',1)
        # 绘制箭头
        data = streamline_sum[j][-1,:]
        x_end,y_end = lonlat_to_xy(data[0],data[1])
        data = streamline_sum[j][-2,:]
        x_endm1,y_endm1 = lonlat_to_xy(data[0],data[1])
        arrow_direction=(np.array([x_end,y_end])-np.array([x_endm1,y_endm1]))
        xy_arrow = np.array([x_end,y_end])
        arrow_color = 'black'
        arrow_width = 0.5
        arrow_0 = np.array([[0,0],[-0.3,0.3],[0.8,0],[-0.3,0.3]])
        # 归一化方向
        a_dn = arrow_direction[:]/xy_ratio[:]
        a_dn = a_dn/(sum(a_dn**2))**0.5
        d = (xy_lim[3]-xy_lim[2]+xy_lim[1]-xy_lim[0])/2
        arrow_1 = arrow_0*arrow_width*0.03*d
        arrow_2 = arrow_1*np.array([[a_dn[0],a_dn[1]],[-a_dn[1],a_dn[0]]])
        xy_ratio_n = xy_ratio/(sum(xy_ratio**2))**0.5             # 归一化比例尺
        arrow_3 = arrow_2*xy_ratio_n+xy_arrow
        xy_ratio_n2 = np.ones((3,1),dtype=float)*xy_ratio_n
        ax.fill(arrow_3[:,0],arrow_3[:,1],color=arrow_color, edgecolor='none')
        # plot_arrow(xy_lim,xy_ratio,[x_end,y_end], [0,0,0],0.5,arrow_direction)

    # 添加网格线 - 设置在最底层
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', zorder=2)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # 保存图像
    filename = f"output/FVCOM/surface_temp_stream_{i_time:02d}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=150)
    plt.close()

    print(f"已保存: {filename}")

print("所有时间步处理完成")