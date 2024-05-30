import numpy as np

# 定义矩阵evec
evec = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9],
                 [10, 11, 12],
                 [13, 14, 15]])

# 使用NumPy的qr函数进行QR分解
Q, R = np.linalg.qr(evec)

# 输出Q和R矩阵
print("Q matrix:")
print(Q)
print("\nR matrix:")
print(R)