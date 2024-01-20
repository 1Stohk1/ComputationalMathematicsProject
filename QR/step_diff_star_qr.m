addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(0, 2); %=======================================================
% Insert data_prep(0, 2) to reproduce the plot with 45 columns
tic
[m, n] = size(A);
atIter = 0;
results = zeros(1, n); 
x_star = A\b;
x_comp = zeros(size(x_star));
for k = 1:n
    atIter = atIter + 1;
    [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
    x_comp(1:k,:) = x;
    results(atIter) = norm(x_comp-x_star)/norm(x_star);
end
toc

semilogy(results, 'g-','LineWidth',2);
hold on

[A, b] = data_prep(0, 0); %=======================================================
% Insert data_prep(0, 2) to reproduce the plot with 45 columns
tic
[m, n] = size(A);
atIter = 0;
results = zeros(1, n); 
x_star = A\b;
x_comp = zeros(size(x_star));
for k = 1:n
    atIter = atIter + 1;
    [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
    x_comp(1:k,:) = x;
    results(atIter) = norm(x_comp-x_star)/norm(x_star);
end
toc

semilogy(results, 'r-','LineWidth',2);
hold on
xlim([0 size(A, 0)])
title('Householder QR')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
% legend('$||Rw-Q^Ty||\slash||y||$', '$||Q^T_2y||\slash||y||$', 'fontsize',16, 'Interpreter', 'Latex')
legend('$\hat{X} ^{(1491\times 45)}$', '$\hat{X} ^{(1491\times 54)}$', '$||Q^T_2y||\slash||y||$','fontsize',16, 'Interpreter', 'Latex')




