addpath(genpath([fileparts(pwd), filesep]));

[A, b] = system_generator(100, 100, 1);

[m, n] = size(A);
atIter = 0;
results = zeros(1, n); 
residuals = zeros(1, n);
times = zeros(1, n);
for k = 1:n
    atIter = atIter + 1;
    for l = 1:1000
        tic;
        [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
        times(atIter) = times(atIter) + toc;
    end
end

times = times  * 1000;

plot(times, 'r-','LineWidth',2);
xlim([0 n])
title('QR Factorization (Time)')
xlabel('Dimension of X')
ylabel('Ex. Time ($\mu$)', 'Interpreter', 'Latex')