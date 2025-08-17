import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pywt

# 读取数据
dataframe = pd.read_csv('./output/HYCOM/flux.txt', sep=' ', header=None, names=['时间','通量'])
data = dataframe['通量'].values

# 小波分析参数设置
wavelet_type = 'db4'  # 使用Daubechies 4小波
decomposition_level = 5  # 分解层数

# 1. 执行小波分解
coeffs = pywt.wavedec(data, wavelet_type, level=decomposition_level)
cA = coeffs[0]  # 近似系数 (最底层)
cD = coeffs[1:]  # 细节系数 (各层)

# 2. 计算各层能量
energy_total = np.sum(data**2)
energy_levels = [np.sum(c**2) for c in coeffs]
energy_ratios = [e/energy_total for e in energy_levels]

# 3. 重构各层信号 (添加长度修正)
reconstructed_levels = []
for i in range(len(coeffs)):
    coeffs_temp = coeffs.copy()
    for j in range(len(coeffs)):
        if j != i:
            coeffs_temp[j] = np.zeros_like(coeffs_temp[j])
    # 重构并截断到原始长度
    rec_signal = pywt.waverec(coeffs_temp, wavelet_type)
    reconstructed_levels.append(rec_signal[:len(data)])

# 4. 阈值去噪 (添加长度修正)
threshold = np.std(data) * np.sqrt(2*np.log(len(data)))  # 通用阈值
coeffs_thresh = [coeffs[0]]  # 保留近似系数
for i in range(1, len(coeffs)):
    coeffs_thresh.append(pywt.threshold(coeffs[i], threshold, mode='soft'))

denoised_data = pywt.waverec(coeffs_thresh, wavelet_type)
denoised_data = denoised_data[:len(data)]  # 关键修正：截断到原始长度

# 5. 绘制结果
plt.figure(figsize=(15, 12))

# 原始信号
plt.subplot(4, 1, 1)
plt.plot(data, 'b')
plt.title('Original Signal')
plt.grid(True)

# 小波系数
plt.subplot(4, 1, 2)
for i, d in enumerate(cD):
    plt.plot(d, label=f'Level {decomposition_level-i}')
plt.plot(cA, 'r', label='Approximation')
plt.title('Wavelet Coefficients')
plt.legend()
plt.grid(True)

# 能量分布
plt.subplot(4, 1, 3)
levels = [f'Approx L{decomposition_level}'] + [f'Detail L{i}' for i in range(decomposition_level,0,-1)]
plt.bar(levels, energy_ratios)
plt.title('Energy Distribution per Level')
plt.ylabel('Energy Ratio')
plt.grid(True)

# 去噪对比
plt.subplot(4, 1, 4)
plt.plot(data, 'b', alpha=0.5, label='Original')
plt.plot(denoised_data, 'r', label='Denoised')
plt.title('Denoising Result')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# 6. 打印分析结果
print("\n小波分析报告:")
print(f"使用小波: {wavelet_type}")
print(f"分解层数: {decomposition_level}")
print(f"总能量: {energy_total:.4f}")
print("\n各层能量占比:")
for i, ratio in enumerate(energy_ratios):
    print(f"{levels[i]}: {ratio*100:.2f}%")

# 计算SNR前确保长度一致
if len(data) == len(denoised_data):
    # 计算信号能量
    signal_energy = np.sum(data**2)
    # 计算噪声能量（即残差平方和）
    noise_energy = np.sum((data - denoised_data)**2)
    # 避免除以0
    if noise_energy == 0:
        snr = float('inf')
    else:
        # 修正括号匹配问题
        snr = 10 * np.log10(signal_energy / noise_energy)
    print(f"\n去噪阈值: {threshold:.4f}")
    print(f"去噪后信噪比(SNR): {snr:.2f} dB")
else:
    print("\n警告: 原始数据和去噪数据长度不一致，无法计算SNR")