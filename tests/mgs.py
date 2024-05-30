import numpy as np

def modified_gram_schmidt(A):
    m, n = A.shape
    Q = np.zeros((m, n))  # 初始化正交矩阵Q

    for j in range(n):
        # 将A的第j列赋值给Q的第j列
        Q[:, j] = A[:, j]

        for i in range(j):
            # 投影Q的第j列到Q的第i列上，并减去这个投影
            Q[:, j] -= Q[:, i] * np.dot(Q[:, i], Q[:, j])

        # 检查是否正交化失败（向量是否为零向量）
        if np.linalg.norm(Q[:, j]) < 1e-10:
            raise ValueError("Zero vector encountered in Modified Gram-Schmidt process.")

        # 归一化
        Q[:, j] /= np.linalg.norm(Q[:, j])

    return Q


def MGS( A ):
    A1 = A[0]
    n,m = A.shape
    q0 =A[0] / np.linalg.norm(A1)
    q = np.zeros((n,m))
    q[0] = q0
    for i in range(n-1):
    # for i in range(m-1):
        for j in range(i+1, n):
        # for j in range(i+1, m):
            A[j] = A[j] - np.inner(A[j],q[i]) * q[i]
        for k in range(m):
            if A[i+1][k] < 1.0e-6:
                q[i+1] = A[i+1]
                '''矩阵的两行相同时，第二个相同元素行会出现一个极小量。若单位化，则会出现[0,0,0,1,0,0,0,0,0,]，
                本步是遍历第二行的元素，若存在极小量，则不需将此行单位化，直接输出该行，那么在结果中需要自己辨别是否存在零向量量'''
            else:
                q[i+1] = A[i+1] / np.linalg.norm(A[i+1])
    return q

# 示例使用
if __name__ == "__main__":
    # 创建一个5x4的矩阵
    # A = np.array([[1, 2, 3, 4],
    #               [5, 6, 7, 8],
    #               [9, 10, 11, 12],
    #               [13, 14, 15, 16],
    #               [17, 18, 19, 20]])
    A = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9],
                 [10, 11, 12],
                 [13, 14, 15]])
    A=A.transpose()
    
    A = np.array([[3,4,5,6],[3,4,5,6],[6,5,3,7]],dtype=np.float64)

    # 调用MGS正交化函数
    # Q = modified_gram_schmidt(A)
    print("A=\n",A)
    Q = MGS(A)

    # 输出正交化后的矩阵
    print("Orthogonalized matrix Q:")
    print(Q)