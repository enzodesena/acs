function obs = drnd(values, pmf, N)
pmf = pmf(:);
thresholds = cumsum(pmf);
uniform_obs = rand(1, N);
obs = zeros(N, 1);
for n = 1:N
    [~, k] = max(uniform_obs(n) < thresholds);
    obs(n) = values(k);
end
end
