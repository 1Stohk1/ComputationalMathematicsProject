addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(0, 2); %=======================================================
% Insert data_prep(0, 2) to reproduce the plot with 45 columns
tic
[m, n] = size(A);
atIter = 0;
results = zeros(1, n); 
x_star = A\b;
% residuals = zeros(1, n);
for k = 1:n
    atIter = atIter + 1;
    [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
%     results(atIter) = norm(A(:, 1:k)*x -b)/norm(b);
%     [Q1, R1] = custom_HQR(A(:, 1:k));
%     x = R1\(Q1*b);
    results(atIter) = norm(x-x_star)/norm(x_star);
%     Rx = R1*x;
%     Qb = Q1*b;
%     residuals(atIter) = norm(Qb(k+1:m,:))/norm(b);
end
toc
semilogy(results, 'g:', 'LineWidth', 2);
hold on
% semilogy(residuals, 'k--', 'LineWidth', 1.5);


[A, b] = data_prep(0, 0); %=======================================================
% Insert data_prep(0, 2) to reproduce the plot with 45 columns
tic
[m, n] = size(A);
atIter = 0;
results = zeros(1, n); 
x_star = A\b;
% residuals = zeros(1, n);
for k = 1:n
    atIter = atIter + 1;
    [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
%     results(atIter) = norm(A(:, 1:k)*x -b)/norm(b);
%     [Q1, R1] = custom_HQR(A(:, 1:k));
%     x = R1\(Q1*b);
    results(atIter) = norm(x-x_star)/norm(x_star);
%     Rx = R1*x;
%     Qb = Q1*b;
%     residuals(atIter) = norm(Qb(k+1:m,:))/norm(b);
end
toc
semilogy(results, 'r-', 'LineWidth', 2);
hold on
% semilogy(residuals, 'k--', 'LineWidth', 1.5);
xlim([0 size(A, 2)])
title('Householder QR')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$\hat{X} ^{(1491\times 45)}$', '$\hat{X} ^{(1491\times 54)}$','fontsize',16, 'Interpreter', 'Latex')
% legend('$||Rw-Q^Ty||\slash||y||$', '$||Q^T_2y||\slash||y||$', 'fontsize',16, 'Interpreter', 'Latex')
