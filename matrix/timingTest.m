function [v,e]= timingTest(A)
tic
[v,e] = eigs(A,100,'SA');
% 154.327414 ��
toc
e=diag(e);
end