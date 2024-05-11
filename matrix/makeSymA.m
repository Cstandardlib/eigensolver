n = 5; % 定义矩阵的大小
symA = zeros(n, n); % 初始化一个 n×n 的零矩阵

% 填充对角线元素
for i = 1:n
    for j = 1:n
        symA(i,j)=i+j-1;
    end
end

tridiagA=tril(triu(symA,1),1);
tridiagA = [
1 2 0 0 0;
2 3 4 0 0;
0 4 5 6 0;
0 0 6 7 8;
0 0 0 8 9];

m=3;
vecs = ones(n, m);

tvecs = tridiagA \ vecs;