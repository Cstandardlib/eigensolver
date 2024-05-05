n = 9; % 矩阵的大小

% 创建矩阵A，其对角线元素从1到9
A = diag(1:n);

% 创建矩阵B，其对角线元素从9到1
B = diag(9:-1:1);

% 使用eig函数求解广义特征值问题 A*x = eigval * B*x
% eigval是包含特征值的向量，x是包含特征向量的列向量
% 注意MATLAB中的eig函数默认返回的是标准的特征值问题，即A*x = eigval*x
% 对于广义特征值问题，MATLAB提供了eig函数的'ch'选项
[eigvec, eigval] = eig(A, B);

% eigval中的特征值是按照升序排列的，我们可以反转它们
eigval = eigval(end:-1:1);
eigvec = eigvec(:, end:-1:1);

% 输出结果
disp('特征值：');
disp(diag(eigval));

disp('对应的特征向量：');
disp(eigvec);