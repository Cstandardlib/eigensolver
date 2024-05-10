n=1000;
a = zeros(n,n); % ��ʼ������aΪn��n�������
dp = 1; % ����dp��ʾ˫���ȣ���MATLAB��Ĭ��Ϊ˫����

for i = 1:n
    a(i,i) = i + 1.0; % ��MATLAB�У�real�������������ʵ������������ֱ��ȡ����ֵ
    for j = 1:(i-1)
        a(j,i) = 1.0 / (i+j);
        a(i,j) = a(j,i); % MATLAB�о����ǶԳƵģ����ֻ��Ҫ��ֵһ��
    end
end

[Va, Ea] = eigs(a,20,'SM');
Ea=diag(Ea);