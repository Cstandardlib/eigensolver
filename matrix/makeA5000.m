n = 5000; % 定义矩阵的大小
a = zeros(n, n); % 初始化一个 n×n 的零矩阵

% 填充对角线元素
for i = 1:n
    a(i,i) = i + 1.0;
    for j = 1:i-1
        a(j,i) = 1.0 / (i+j);
        a(i,j) = a(j,i);
    end
end

b=diag(ones(n, 1) * 2);

tic
[Va, Ea] = eigs(a,b,1,'SA'); % 92.057098秒
% [Va, Ea] = eigs(a,b,200,'SA'); % 92.057098秒
toc
ea=diag(Ea);