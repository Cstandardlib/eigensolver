function [v,e]= timingTest(A)
tic
[v,e] = eigs(A,199,'SM'); % 154.327414 √Î
toc
end