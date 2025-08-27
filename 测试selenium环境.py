import undetected_chromedriver as uc

options = uc.ChromeOptions()
# 指定 Chromium 浏览器的可执行文件路径
options.binary_location = "/usr/bin/chromium"  # 请根据你的系统实际情况修改路径

# 可以同时指定用户数据目录（可选）
# options.user_data_dir = "/path/to/your/chromium/profile"

try:
    # 初始化 driver，注意新版本可能不再使用 executable_path 参数 :cite[8]
    driver = uc.Chrome(options=options)
    options.add_argument("--remote-debugging-port=9222")

    # 打开网页进行测试
    driver.get("https://www.baidu.com")
    print(driver.page_source)  # 打印页面内容，确认成功访问
    driver.quit()
    print("测试成功！")
except Exception as e:
    print(f"发生错误: {e}")