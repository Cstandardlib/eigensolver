n = 9; % ����Ĵ�С

% ��������A����Խ���Ԫ�ش�1��9
A = diag(1:n);

% ��������B����Խ���Ԫ�ش�9��1
B = diag(9:-1:1);

% ʹ��eig��������������ֵ���� A*x = eigval * B*x
% eigval�ǰ�������ֵ��������x�ǰ�������������������
% ע��MATLAB�е�eig����Ĭ�Ϸ��ص��Ǳ�׼������ֵ���⣬��A*x = eigval*x
% ���ڹ�������ֵ���⣬MATLAB�ṩ��eig������'ch'ѡ��
[eigvec, eigval] = eig(A, B);

% eigval�е�����ֵ�ǰ����������еģ����ǿ��Է�ת����
eigval = eigval(end:-1:1);
eigvec = eigvec(:, end:-1:1);

% ������
disp('����ֵ��');
disp(diag(eigval));

disp('��Ӧ������������');
disp(eigvec);