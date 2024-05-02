function ea = timingTest(A)
tic
ea = eigs(A,199,'SM'); % 154.327414 √Î
toc
end