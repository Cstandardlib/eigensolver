import matplotlib.pyplot as plt

# x axis: different matrices
# iterations = [0, 1, 2, 3]
iterations = ["Si2","Na5","Si5H12","Ga3As3H12"]


# 每组数据的矩阵值，这里只是示例数据，你需要替换成你自己的数据
# iter for ["Si2","Na5","Si5H12","Ga3As3H12"]
# thinQR,         nopre   diag   tri
# dense5000,        80     16   
# Si2 200-225       12     12     11
# Na5 100-120       49     50     40
# Na5 200-225       45     46     38
# Si5H12 200-225    74     73     60
# Ga3As3H12 10-20   854    658    388
nopre_matrix_values = [0, 10, 20, 25]
diag_matrix_values = [1, 15, 28, 35]
tri_matrix_values = [2, 12, 22, 30]

# 设置图表
plt.figure(figsize=(10, 6))

# 绘制每组数据的曲线
plt.plot(iterations, nopre_matrix_values, marker='o', label='No Preconditioner')
plt.plot(iterations, diag_matrix_values, marker='s', label='Diagonal Preconditioner')
plt.plot(iterations, tri_matrix_values, marker='^', label='Triangular Preconditioner')

# 添加图例
plt.legend()

# 添加标题和轴标签
plt.title('Iter-Matrix Curves for Different Conditions')
plt.xlabel('Iteration')
plt.ylabel('Matrix Value')

# 显示网格
plt.grid(True)

# 显示图表
plt.show()