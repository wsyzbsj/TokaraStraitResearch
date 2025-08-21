import numpy as np

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
    xy_start = np.interp2(xn,yn,mask,zxx,zyy)
    in = find(xy_start>0)
    xy_start(in) = 1
    zy = np.linspace(0,1,xy_end.shape[0])
    zx = np.linspace(0,1,xy_end.shape[1])
    [zxx,zyy] = np.meshgrid(zx,zy)
    xy_end = np.interp2(xn,yn,mask,zxx,zyy)
    in = find(xy_end>0)
    xy_end(in) = 1

    # 将流线划分为num种，速度越大的流线越长
    length_sl=linspace(5,40,num)
    V2=(un**2+vn**2)**0.5
    V2_max=V2.max()
    V2_min=V2.min()

    V2_space=np.linspace(V2_min,V2_max,num+1)

    # 1当xy_start内还有可放置的新点的位置时，进行循环
    k=0 #循环次数,也是流线个数
    while ~all(xy_start,'all'):
        k=k+1;
        # 2随机一个start内网格点作为种子点
        [start_id_y,start_id_x]=find(xy_start==0)
        randnum=randi(start_id_y.shape[0])
        x_pos_i=id2axis(dstart,start_id_x(randnum,1))
        y_pos_i=id2axis(dstart,start_id_y(randnum,1))
        streamline_seed(k,:)=[x_pos_i,y_pos_i]                                                      # 保存种子点
        V2_seed=interp2(xn,yn,V2,x_pos_i,y_pos_i)                                                   # 计算种子点处的速度
        [~,~,sl_N] = histcounts(V2_seed,V2_space)
        if sl_N==0 || isnan(sl_N)
            sl_N =5;
        end

        num_streamline=round(length_sl(sl_N)) * magnify;

        # 3绘制流线
        streamline_i_1 = stream2(xn,yn, un, vn,x_pos_i,y_pos_i,[0.1,num_streamline])
        streamline_i_2 = stream2(xn,yn,-un,-vn,x_pos_i,y_pos_i,[0.1,num_streamline])
        # 4以xy_end为标准，删除自相交或间隔太近的点。并顺便标记xy_end
        [streamline_i_1,xy_end,xy_start]=delete_self(streamline_i_1{1}, xy_end,dend,xy_start,dstart)
        [streamline_i_2,xy_end,xy_start]=delete_self(streamline_i_2{1}, xy_end,dend,xy_start,dstart)
        # 5保存
        streamline_k=[flipud(streamline_i_2);streamline_i_1(2:end,:)]                               # 新的流线
        streamline_sum{k}=[xmin+streamline_k(:,1)*(xmax-xmin), ymin+streamline_k(:,2)*(ymax-ymin)]  # 从归一化还原
    end
    streamline_seed=[streamline_seed(:,1)*(xmax-xmin)+xmin, streamline_seed(:,2)*(ymax-ymin)+ymin]
    end



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