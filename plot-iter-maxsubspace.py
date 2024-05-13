import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# 假设文件名为 'lobpcg_output.txt'
filename = 'Ga3As3H12cho_1.txt'

# Si2, to seek 10 pairs with different n_max_subspace

# 初始化变量
# Si2 10-subspace
# n_max_subspace_list =  list(range(15, 51, 5)) #[i for i in ]
# iter_list = [38,30,27,24,22,20,19,18]
# Na5
n_max_subspace_list =  list(range(15, 51, 5)) #[i for i in ]
iter_list = [38,30,27,24,22,20,19,18]


# # 读取文件并解析每一行
# with open(filename, 'r') as file:
#     for line in file:
#         # 使用正则表达式匹配一行中的五个数值
#         match = re.match(r'([-+]?\d+)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)', line)
#         if match:
#             n_max_subspace_list.append(float(match.group(1)))  # 提取第一个数
#             iter_list.append(float(match.group(4)))  # 提取第四个数

# n_max_subspace_list = n_max_subspace_list[0:60]
# iter_list = iter_list[0:60]
# 打印提取的数值以验证
print("First number list:", n_max_subspace_list)
print("Fourth number list:", iter_list)

# 如果需要绘图或其他处理，可以继续使用这两个列表
# 绘图
plt.figure(figsize=(10, 5))
plt.plot(n_max_subspace_list, iter_list, marker='o')
plt.title('LOBPCG Iterations vs n_space')
plt.xlabel('n_space')
# plt.title('LOBPCG Iterations vs n_max_subspace')
# plt.xlabel('n_max_subspace')
plt.ylabel('Iterations')
plt.grid(True)

# 设置X轴和Y轴的刻度为整数
ax = plt.gca()  # gca代表get current axis，获取当前的Axes对象
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))

plt.show()