import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


# 初始化变量
# Si2 10-subspace
# n_max_subspace_list =  list(range(15, 51, 5)) #[i for i in ]
# iter_list = [38,30,27,24,22,20,19,18]
# Si2
n_max_subspace_list =  list(range(30, 37)) #[i for i in ]
iter_list_no = [77,75,60,49,48,43,41]
iter_list_diag = [77,64,53,49,43,42,39]
iter_list_tri = [71,50,48,42,37,35,33]
# Na5
# n_max_subspace_list =  list(range(100, 108)) #[i for i in ]
# iter_list_no = [149,130,104,92,85,81,76,71]
# iter_list_diag = [150,127,95,92,82,80,73,71]
# iter_list_tri = [130,96,87,83,71,66,63,58]

# 如果需要绘图或其他处理，可以继续使用这两个列表
# 绘图
plt.figure(figsize=(10, 5))
plt.figure(figsize=(10, 5))
plt.plot(n_max_subspace_list, iter_list_no, marker='o', label='No Preconditioner', color='blue')
plt.plot(n_max_subspace_list, iter_list_diag, marker='o', label='Diagonal', color='red')#, linestyle='--')
plt.plot(n_max_subspace_list, iter_list_tri, marker='o', label='Tridiagonal', color='green')#, linestyle=':')
# plt.plot(n_max_subspace_list, iter_list_no, marker='o', label='No Preconditioner')
# plt.plot(n_max_subspace_list, iter_list_diag, marker='o', label='Diagonal')
# plt.plot(n_max_subspace_list, iter_list_tri, marker='o', label='Tridiagonal')

# 添加图例
plt.legend()

plt.title('Si2 LOBPCG Iterations vs n_space')
plt.xlabel('n_space')
plt.ylabel('Iterations')
plt.grid(True)

# 设置X轴和Y轴的刻度为整数
ax = plt.gca()  # gca代表get current axis，获取当前的Axes对象
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))

plt.show()