

%% Testing drnd(values, pmf, N)
N = 100000;
values = [-1, -2, -3];
pmf = [0.1, 0.3, 0.6];
obs = drnd(values, pmf, N);
for i=1:length(values)
    epmf(i) = sum(obs==values(i))/N;
end

assert(all(abs(epmf-pmf) < 0.01));

display('All tests passed!');