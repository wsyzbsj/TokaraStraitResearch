# /usr/bin/python
coding = 'utf-8'

import requests
from bs4 import UnicodeDammit
from lxml import html
import os,pathlib,glob,openpyxl

string_html_url_prefix = "https://www.jodc.go.jp/TideAttrib/tide_attrib_"
string_html_url_subfix = "_e.htm"

# 输出表格
workbook = openpyxl.Workbook()

# 读取站号
storage_paths = sorted(glob.glob(r'D:\Documents\Data\JODC\验潮站数据\*'))
for i_storage_path,storage_path in enumerate(storage_paths):
    if os.path.isdir(storage_path):
        sta_id = pathlib.Path(storage_path).stem
        print(sta_id)
        if i_storage_path == 0:
            worksheet = workbook.active
        else:
            worksheet = workbook.create_sheet()
        worksheet.title = sta_id
        url = string_html_url_prefix + sta_id + string_html_url_subfix
        string_resource = requests.get(url)
        string_html_content = string_resource.content
        dammit = UnicodeDammit(string_html_content,
                               ["utf-16-le", "shift-jis", "utf-8", "euc-jp", "cp932", "iso-8859-1"],  # 优先尝试的编码列表
                               is_html=True)
        decoded_content = dammit.unicode_markup
        # 转换为XPath可解析对象
        tree = html.fromstring(decoded_content)
        # 使用XPath获取所有td
        tds = tree.xpath('//tr/td | //tr/TD | //TR/td | //TR/TD')  # 处理大小写问题
        # 按行分组
        rows = {}
        for td in tds:
            tr = td.getparent()
            if tr not in rows:
                rows[tr] = []
            rows[tr].append(td)
        # 遍历处理
        for i_tr, (tr, tds) in enumerate(rows.items()):
            for i_td, td in enumerate(tds):
                try:
                    name = td.text_content().replace('\r\n','')
                except AttributeError:
                    name = ''
                # 预期结果：
                #| Station Code | xxx |
                #| Station Name | xxx |
                # 空行
                #| Date | Career |
                #| xxxx | xxxxxx |
                # 处理缺失
                if len(tds) ==1:
                    worksheet.cell(row=i_tr+1,column=2).value = name
                else:
                    worksheet.cell(row=i_tr+1,column=i_td+1).value = name
    else:
        continue

workbook.save('验潮站标准水位.xlsx')