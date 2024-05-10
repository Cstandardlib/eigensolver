function [v,e]= timingTest(A)
tic
[v,e] = eigs(A,1,'SM'); % 154.327414 √Î
toc
end