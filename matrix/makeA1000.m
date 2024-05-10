n=1000;
a = zeros(n,n); % 初始化矩阵a为n×n的零矩阵
dp = 1; % 假设dp表示双精度，在MATLAB中默认为双精度

for i = 1:n
    a(i,i) = i + 1.0; % 在MATLAB中，real函数返回输入的实部，这里用来直接取整数值
    for j = 1:(i-1)
        a(j,i) = 1.0 / (i+j);
        a(i,j) = a(j,i); % MATLAB中矩阵是对称的，因此只需要赋值一次
    end
end

[Va, Ea] = eigs(a,20,'SM');
Ea=diag(Ea);