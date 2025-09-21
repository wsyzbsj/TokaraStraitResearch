import numpy as np
import pandas as pd
import os
import glob

storage_paths = '/home/yzbsj/Data/海洋数据/JODC/验潮站数据'

# 读取所有Sheet
dfs = list(pd.read_excel('不同站点修正数据.xlsx', sheet_name=None))
for df in dfs:
    dataframe_fix = pd.read_excel('不同站点修正数据.xlsx', sheet_name=df)
    column_names = dataframe_fix.columns.to_list()
    column_names.remove('Date')
    column_names.remove('Latitude')
    column_names.remove('Longitude')

    # 读取水位信息
    storage_path =os.path.join(storage_paths,df)
    files = sorted(glob.glob(os.path.join(storage_path,'*.txt')))
    with open('output/验潮站数据/'+df+'.csv','w') as fout:
        fout.write('站号,时间,水位\n')
        for file in files:
            print(file)
            with open(file,'r') as fin:
                # 获取总长度
                text = fin.read()
                char_count = len(text)
                fin.seek(0)
                i_count = 0
                i_hour = 0
                date = ''
                data = ''
                while i_count < char_count:
                    temp = fin.read(1)
                    if temp == ',':
                        if i_hour == 0: # 站号
                            sta_id = data
                        elif i_hour ==1:# 日期
                            date = data
                        else:   # 0、1为站号和日期
                            fout.write(sta_id+','+date+' '+f"{i_hour-2:02d}"+":00:00"+','+data+'\n')
                        data = ''
                        i_hour += 1     # 小时+1(-2为最终真实时间)
                        i_count +=1
                    elif temp == '\n':
                        fout.write(sta_id+','+date+' '+f"{i_hour-2:02d}"+":00:00"+','+data+'\n')
                        data = ''
                        date = ''
                        i_hour = 0      # 换行表示时间为0时
                        i_count +=1
                    else:
                        data += temp
                        i_count +=1
        fout.close()

    # 修正
    dataframe = pd.read_csv('output/验潮站数据/'+df+'.csv')
    dataframe = dataframe.drop(dataframe[dataframe['水位'] == 9999].index)        # 删除无效数据
    dataframe = dataframe.reset_index(drop=True)
    dates_fix = dataframe_fix['Date'].to_numpy()                                # 获取修正时间
    dates = dataframe_fix['时间'].to_numpy()                                      # 获取验潮站观测时间
    # 计算位置
    loc = []
    loc.append(0)
    for date_fix in dates_fix:
        loc.append(np.searchsorted(dates, date_fix))
    for i in range(len(loc)):
        begin = loc
        if i==len(loc)-1:
            end = len(dates)
        else:
            end = loc[i+1]
        for j in range(begin, end):
            # 根据时间修正水位数据
            i_time = 0
            date = dataframe_fix.loc[j,'时间']
            for column_name in column_names:
                dataframe.loc[j,'水位'+'_'+column_name] = dataframe.loc[j,'水位']+dataframe_fix.loc[i,column_name]
    dataframe.to_csv('output/验潮站数据/'+df+'_fix.csv')