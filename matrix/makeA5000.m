n = 5000; % �������Ĵ�С
a = zeros(n, n); % ��ʼ��һ�� n��n �������

% ���Խ���Ԫ��
for i = 1:n
    a(i,i) = i + 1.0;
    for j = 1:i-1
        a(j,i) = 1.0 / (i+j);
        a(i,j) = a(j,i);
    end
end

b=diag(ones(n, 1) * 2);

tic
[Va, Ea] = eigs(a,b,1,'SA'); % 92.057098��
% [Va, Ea] = eigs(a,b,200,'SA'); % 92.057098��
toc
ea=diag(Ea);