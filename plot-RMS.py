import re
import matplotlib.pyplot as plt

# 假设文件名为 'lobpcg_output.txt'
# filename = 'Ga3As3H12cho_1.txt'
# filename = 'out/gen1000thin_1_3_nopre.txt'
# filename = 'out/Si2_1_pre_cho.txt'
# filename = 'out/Si2thin_1_2_diag.txt'
# filename = 'out/Na5cho_1_nopre.txt'
# filename = 'out/Na5cho_1.txt'
# filename = 'out/Si5H12cho_1.txt'
# filename = 'out/Ga3As3H12cho_1.txt'
# filename = 'data/Ga3As3H12thin_1_2_nopre.txt'
# filename = 'data/Ga3As3H12thin_1_2_diag.txt'
filename = 'data/Ga3As3H12thin_1_2_tri.txt'


# 初始化变量
first_num_list = []
fourth_num_list = []

# 读取文件并解析每一行
with open(filename, 'r') as file:
    for line in file:
        # 使用正则表达式匹配一行中的五个数值
        match = re.match(r'([-+]?\d+)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)', line)
        if match:
            first_num_list.append(float(match.group(1)))  # 提取第一个数
            fourth_num_list.append(float(match.group(4)))  # 提取第四个数

# first_num_list = first_num_list[0:30]
# fourth_num_list = fourth_num_list[0:30]
# Ga
first_num_list = first_num_list[0:60]
fourth_num_list = fourth_num_list[0:60]
# 打印提取的数值以验证
print("First number list:", first_num_list)
print("Fourth number list:", fourth_num_list)

# 如果需要绘图或其他处理，可以继续使用这两个列表
# 绘图
plt.figure(figsize=(10, 5))
plt.plot(first_num_list, fourth_num_list, marker='o')
plt.title('LOBPCG Residuals vs Iterations')
plt.xlabel('Iteration')
plt.ylabel('Residuals')
plt.grid(True)
plt.show()