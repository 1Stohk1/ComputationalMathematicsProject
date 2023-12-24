addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(0, 0); %=======================================================
tic
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
%     results(atIter) = norm(A(:, 1:k)*x -b)/norm(b);
%     [Q1, R1] = custom_HQR(A(:, 1:k));
%     results(atIter) = norm(R1*x-Q1*b)/norm(b);
%     Rx = R1*x;
%     Qb = Q1*b;
%     residuals(atIter) = norm(Qb(k+1:m,:))/norm(b);
end
toc
plot(results, 'LineWidth', 2);
hold on
plot(residuals, 'r--', 'LineWidth', 1.5);

title('Householder QR')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$||Rw-Q^Ty||\slash||y||$', '$||Q^T_2y||\slash||y||$', 'fontsize',16, 'Interpreter', 'Latex')

clf;

times = times  * 1000;

plot(times, 'r-','LineWidth',2);
xlim([0 55])
title('QR Factorization (Time)')
xlabel('Iterations')
ylabel('Ex. Time ($\mu$)', 'Interpreter', 'Latex')