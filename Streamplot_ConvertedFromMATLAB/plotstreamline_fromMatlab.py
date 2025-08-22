# 注:从matlab翻译,stream2更改为from scipy.integrate import odeint
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import RegularGridInterpolator

def id2axis(distance,id):
    n=np.ceil((0.5-distance/2)/distance)*2+1        # 分割的数量
    min_distance=(1-(n-2)*distance)/2               # 两端最小的距离
    if id==1:
        xpoint=min_distance/2
    elif id==n:
        xpoint=1-min_distance/2
    else:
        xpoint=min_distance+(id-1.5)*distance
    return xpoint

    # 原MATLAB代码:
    # function xpoint=id2axis(distance,id)
    # %取网格的中点
    #
    # N=ceil((0.5-distance/2)/distance)*2+1;%分割的数量
    # min_distance=(1-(N-2)*distance)/2;%两端最小的距离
    # if id==1
    #     xpoint=min_distance/2;
    # elseif id==N
    # xpoint=1-min_distance/2;
    # else
    # xpoint=min_distance+(id-1.5)*distance;
    # end
    # end

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
    N=np.ceil((0.5-distance/2)/distance)*2+1                                                                            # 分割的数量
    min_distance=(1-(N-2)*distance)/2                                                                                   # 两端最小的距离

    # x的位置
    if x<=min_distance:
        pos_id_x=1
    elif x>=1-min_distance:
        pos_id_x=N
    else:
        pos_id_x=np.ceil((x-min_distance)/distance)+1

    # y的位置
    if y<=min_distance:
        pos_id_y=1
    elif y>=1-min_distance:
        pos_id_y=N
    else:
        pos_id_y=np.ceil((y-min_distance)/distance)+1

    # xy转ind
    pos_id=sub2ind([N,N],pos_id_y,pos_id_x)

    return pos_id

    # 原MATLAB代码:
    # function pos_id=axis2id(x,y,distance)
    #
    # N=ceil((0.5-distance/2)/distance)*2+1;%分割的数量
    # min_distance=(1-(N-2)*distance)/2;%两端最小的距离
    #
    # %x的位置
    # if x<=min_distance
    #     pos_id_x=1;
    # elseif x>=1-min_distance
    # pos_id_x=N;
    # else
    # pos_id_x=ceil((x-min_distance)/distance)+1;
    # end
    #
    # %y的位置
    # if y<=min_distance
    #     pos_id_y=1;
    # elseif y>=1-min_distance
    # pos_id_y=N;
    # else
    # pos_id_y=ceil((y-min_distance)/distance)+1;
    # end
    #
    # %xy转ind
    # pos_id=sub2ind([N,N],pos_id_y,pos_id_x);
    # end

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

    # 原MATLAB代码:
    # function [sl_i,xy_end,xy_start]=delete_self(sl_i,xy_end,dend,...
    # xy_start,dstart)
    # %lidd rm NaN
    # in = isnan(sl_i(:,1));
    # sl_i = sl_i(~in,:);
    #
    # %sl_i streamline流线，两列N行形式
    # N=size(sl_i,1);
    # pos_id_last=axis2id(sl_i(1,1),sl_i(1,2),dend);
    # xy_end(pos_id_last)=1;%第一个点标记
    #
    # %顺便标记xy_start
    # pos_id_s=axis2id(sl_i(1,1),sl_i(1,2),dstart);
    # xy_start(pos_id_s)=1;
    # for j=2:N
    # pos_id_now=axis2id(sl_i(j,1),sl_i(j,2),dend);
    # if pos_id_now~=pos_id_last
    # %如果现在的点和原有的点在同一区域，则不管它
    # %如果不在同一区域，检测新的点是否已经被占用
    # if xy_end(pos_id_now)==1
    #     %如果该点被占用，说明出现与其它流线太近的情况，则直接停止
    #     j=j-1;
    #     break
    # else
    #     %如果没被占用，则把新点添加上
    #     xy_end(pos_id_now)=1;
    #     pos_id_last=pos_id_now;
    # end
    # end
    # %顺便标记xy_start
    # pos_id_s=axis2id(sl_i(j,1),sl_i(j,2),dstart);
    # xy_start(pos_id_s)=1;
    # end
    # sl_i(j:end,:)=[];
    # end

# 模拟 MATLAB 的 histcounts 函数，只获取 bin 索引
def matlab_histcounts_indices(data, bins):
    # 使用 digitize 获取 bin 索引
    indices = np.digitize(data, bins)

    # 处理超出范围的值（MATLAB 中返回 0）
    # 对于小于最小 bin 的值，digitize 返回 0，与 MATLAB 一致
    # 对于大于最大 bin 的值，digitize 返回 len(bins)+1，需要调整为 0
    indices[indices > len(bins)] = 0

    return indices

def plot_arrow(xy_lim,xy_ratio,xy_arrow,arrow_color,arrow_width,arrow_direction) -> None:
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
    plt.fill(arrow_3[:,0],arrow_3[:,1],color=arrow_color, edgecolor='none')
    # 下附原MATLAB代码:
    #
    # function plot_arrow(xy_lim,xy_ratio,xy_arrow,arrow_color,arrow_width, ...
    # arrow_direction)
    #
    # %初始化箭头形状（归一化的形状）
    # % arrow_0=[0,0;-1,0.5;-1,-0.5];
    # % arrow_0=[1, 0  ;
    # %          0, 0.5;
    # %          0,-0.5];
    # arrow_0=[ 0   , 0  ;
    # -0.3, 0.3;
    # 0.8   , 0  ;
    # -0.3,-0.3];
    # %对方向进行归一化
    # a_dn=arrow_direction(:)./xy_ratio(:);
    # a_dn=a_dn/sqrt(sum(a_dn.^2));
    # d=(xy_lim(4)-xy_lim(3)+xy_lim(2)-xy_lim(1))/2;
    # %箭头对窗口缩放
    # arrow_1=arrow_0*arrow_width*0.03*d;
    # %箭头旋转
    # arrow_2=arrow_1*[a_dn(1),a_dn(2);-a_dn(2),a_dn(1)];
    # %箭头变形
    # xy_ratio_n=xy_ratio/sqrt(sum(xy_ratio.^2));%对比例尺归一化
    # %arrow_3=arrow_2.*xy_ratio_n+xy_arrow;%由于兼容性问题，更改为下面语句
    # % xy_ratio_n2=ones(3,1)*xy_ratio_n;
    # xy_ratio_n2=ones([size(arrow_0,1),1])*xy_ratio_n;
    # arrow_3=arrow_2.*xy_ratio_n2+xy_arrow;
    # fill(arrow_3(:,1),arrow_3(:,2),arrow_color,'EdgeColor','none')
    # end

def my_streamline_mutli(x,y,u,v,dstart,num,magnify) -> list:
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

    num_start=np.ceil((0.5-dstart/2)/dstart)*2+1
    num_end=np.ceil((0.5-dend/2)/dend)*2+1

    # 初始化所有网格点，0代表可以放置新点，1代表已经存在原有的点
    xy_start = np.zeros((num_start,num_start),dtype=float)
    xy_end = np.zeros((num_end,num_end),dtype=float)

    # 标记陆地上点为1:，代表不会进入循环
    mask = np.isnan(un).astype(float)
    zy = np.linspace(0,1,xy_start.shape[0])
    zx = np.linspace(0,1,xy_start.shape[1])
    [zxx,zyy] = np.meshgrid(zx,zy)
    # 原本插值:xy_start = np.interp2(xn,yn,mask,zxx,zyy)
    # 创建插值器
    interpolator = RegularGridInterpolator((yn[:, 0], xn[0, :]), mask, method='linear', bounds_error=False, fill_value=np.nan)
    # 进行插值
    xy_start = interpolator((zyy, zxx))  # 注意坐标顺序

    py_in = np.where(xy_start.flatten(order='C')>0)
    xy_start[py_in] = 1
    zy = np.linspace(0,1,xy_end.shape[0])
    zx = np.linspace(0,1,xy_end.shape[1])
    [zxx,zyy] = np.meshgrid(zx,zy)
    # 原本的插值:xy_end = np.interp2(xn,yn,mask,zxx,zyy)
    # 创建插值器
    interpolator = RegularGridInterpolator((yn[:, 0], xn[0, :]), mask, method='linear', bounds_error=False, fill_value=np.nan)
    # 进行插值
    xy_end = interpolator((zyy, zxx))  # 注意坐标顺序
    py_in = np.where(xy_end.flatten(order='C')>0)
    xy_end[py_in] = 1

    # 将流线划分为num种，速度越大的流线越长
    length_sl=np.linspace(5,40,num)                                                                                     # 按速度分类
    V2=(un**2+vn**2)**0.5
    V2_max=V2.max()
    V2_min=V2.min()

    V2_space=np.linspace(V2_min,V2_max,num+1)

    # 1. 当xy_start内还有可放置的新点的位置时，进行循环
    k=0                                                                                                                 # 流线数(循环次数)
    streamline_seed = []
    while not np.all(xy_start):
        # 2. 随机一个start内网格点作为种子点
        [start_id_y,start_id_x]=np.where(xy_start==0)
        randnum=np.random.randint(1, start_id_y.shape[0]+1)
        x_pos_i=id2axis(dstart,start_id_x[randnum,0])
        y_pos_i=id2axis(dstart,start_id_y[randnum,0])
        streamline_seed[k].append([x_pos_i,y_pos_i])                                                                    # 保存种子点
        # 原本的插值:V2_seed=np.interp2(xn,yn,V2,x_pos_i,y_pos_i)                                                            # 计算种子点处的速度
        # 创建插值器
        interpolator = RegularGridInterpolator((yn[:, 0], xn[0, :]), V2, method='linear', bounds_error=False, fill_value=np.nan)
        # 进行插值
        V2_seed = interpolator((zyy, zxx))  # 注意坐标顺序
        sl_N = matlab_histcounts_indices(V2_seed, V2_space)
        # [~,~,sl_N] = histcounts(V2_seed,V2_space)
        if sl_N==0 or np.isnan(sl_N):
            sl_N =5

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
        interpolator_u = RegularGridInterpolator((yn[:,0], xn[0,:]), un.T)
        interpolator_v = RegularGridInterpolator((yn[:,0], xn[0,:]), vn.T)

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
        streamline_sum[k]=[xmin+streamline_k[:][0]*(xmax-xmin), ymin+streamline_k[:][1]*(ymax-ymin)]                      # 从归一化还原
        k += 1
    streamline_seed=[streamline_seed[:][0]*(xmax-xmin)+xmin, streamline_seed[:][1]*(ymax-ymin)+ymin]

    return [streamline_sum,streamline_seed]

    # 下附MATLAB原代码:
    # function [streamline_sum,streamline_seed]=my_streamline_mutli(x,y,u,v, ...
    # dstart,num,magnify)
    # %0处理前设置
    # %设置网格密度(01区间内归一化的长度）
    # %dstart=0.05;默认0.05
    # dend=0.5*dstart;
    #
    # %xmin=min(x,[],'all');xmax=max(x,[],'all');
    # %ymin=min(y,[],'all');ymax=max(y,[],'all');
    # xmin=min(min(min(x)));xmax=max(max(max(x)));
    # ymin=min(min(min(y)));ymax=max(max(max(y)));
    #
    # %归一化，将流场缩放为0-1区间的矩形
    # xn=(x-xmin)/(xmax-xmin);
    # yn=(y-ymin)/(ymax-ymin);
    # un=u/(xmax-xmin);
    # vn=v/(ymax-ymin);
    #
    # num_start=ceil((0.5-dstart/2)/dstart)*2+1;
    # num_end=ceil((0.5-dend/2)/dend)*2+1;
    #
    # %初始化所有网格点，0代表可以放置新点，1代表已经存在原有的点
    # xy_start=zeros(num_start,num_start);
    # xy_end=zeros(num_end,num_end);
    #
    # % 标记陆地上点为1:，代表不会进入循环
    # mask = double(isnan(un));
    # zy = linspace(0,1,size(xy_start,1));
    # zx = linspace(0,1,size(xy_start,2));
    # [zxx,zyy] = meshgrid(zx,zy);
    # xy_start = interp2(xn,yn,mask,zxx,zyy);
    # in = find(xy_start>0);
    # xy_start(in) = 1;
    # zy = linspace(0,1,size(xy_end,1));
    # zx = linspace(0,1,size(xy_end,2));
    # [zxx,zyy] = meshgrid(zx,zy);
    # xy_end = interp2(xn,yn,mask,zxx,zyy);
    # in = find(xy_end>0);
    # xy_end(in) = 1;
    #
    # %将流线划分为num种，速度越大的流线越长
    # length_sl=linspace(5,40,num);
    # V2=sqrt(un.^2+vn.^2);
    # V2_max=max(max(max(V2)));
    # V2_min=min(min(min(V2)));
    # %V2_min=min(V2,[],'all');V2_max=max(V2,[],'all');
    #
    # V2_space=linspace(V2_min,V2_max,num+1);
    #
    # %1当xy_start内还有可放置的新点的位置时，进行循环
    # k=0;%循环次数，也是流线个数
    # while ~all(xy_start,'all')
    #     k=k+1;
    #     %2随机一个start内网格点作为种子点
    #     [start_id_y,start_id_x]=find(xy_start==0);
    #     randnum=randi(size(start_id_y,1));
    #     x_pos_i=id2axis(dstart,start_id_x(randnum,1));
    #     y_pos_i=id2axis(dstart,start_id_y(randnum,1));
    #     streamline_seed(k,:)=[x_pos_i,y_pos_i];%保存种子点
    #     V2_seed=interp2(xn,yn,V2,x_pos_i,y_pos_i);%计算种子点处的速度
    #     [~,~,sl_N] = histcounts(V2_seed,V2_space);
    #     if sl_N==0 || isnan(sl_N)
    #         sl_N =5;
    #     end
    #
    #     num_streamline=round(length_sl(sl_N)) * magnify;
    #
    #     %3绘制流线
    #     streamline_i_1 = stream2(xn,yn, un, vn,x_pos_i,y_pos_i,[0.1,num_streamline]);
    #     streamline_i_2 = stream2(xn,yn,-un,-vn,x_pos_i,y_pos_i,[0.1,num_streamline]);
    #     %4以xy_end为标准，删除自相交或间隔太近的点。并顺便标记xy_end
    #     [streamline_i_1,xy_end,xy_start]=delete_self(streamline_i_1{1}, ...
    #     xy_end,dend,xy_start,dstart);
    #     [streamline_i_2,xy_end,xy_start]=delete_self(streamline_i_2{1}, ...
    #     xy_end,dend,xy_start,dstart);
    #     %5保存
    #     streamline_k=[flipud(streamline_i_2);streamline_i_1(2:end,:)];%新的流线
    #     streamline_sum{k}=[xmin+streamline_k(:,1)*(xmax-xmin), ...
    #     ymin+streamline_k(:,2)*(ymax-ymin)];%从归一化还原
    # end
    # streamline_seed=[streamline_seed(:,1)*(xmax-xmin)+xmin, ...
    # streamline_seed(:,2)*(ymax-ymin)+ymin];
    # end