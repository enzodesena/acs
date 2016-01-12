function run_tests

%% Testing drnd(values, pmf, N)
N = 100000;
values = [-1, -2, -3];
pmf = [0.1, 0.3, 0.6];
obs = drnd(values, pmf, N);
for i=1:length(values)
    epmf(i) = sum(obs==values(i))/N;
end
assert(is_equal(epmf, pmf, 0.01));

%% Testing vMUM standard
[mu, k, p1, p2, p3] = circ_vmum_standard(pi, 10, 0.3, 0.5, 0.2);
assert(is_equal(circ_vmum_pdf(mu, k, p1, p2, p3, linspace(0,2*pi,100)),...
                circ_vmum_pdf(pi, 10, 0.3, 0.5, 0.2, linspace(0,2*pi,100))));
            
[mu, k, p1, p2, p3] = circ_vmum_standard(pi, -10, 0.3, 0.5, 0.2);
assert(is_equal(circ_vmum_pdf(mu, k, p1, p2, p3, linspace(0,2*pi,100)),...
                circ_vmum_pdf(pi, -10, 0.3, 0.5, 0.2, linspace(0,2*pi,100))));

display('All tests passed!');

end

function result = is_equal(a, b, tol)
if nargin == 2
    tol = 1e-10;
end

assert(all(abs(epmf-pmf) < 0.01));

%% Testing vM for large kappa
thetas = linspace(-4*pi,4*pi,1000);
vm_uncorrected = circ_vm_pdf(0, 600, thetas, false);
vm_corrected = circ_vm_pdf(0, 600, thetas, true);
assert(all(abs(vm_corrected-vm_uncorrected)<1e-2));

assert(all(not(isnan(circ_vm_pdf(0, 10000, thetas, true)))));

display('All tests passed!');