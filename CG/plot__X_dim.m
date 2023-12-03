addpath(genpath([fileparts(pwd), filesep]));


[A, b] = data_prep(1, -1); %=======================================================
[x, res] = custom_conjgrad(A, b, b);
result = norm(A*x -b)/norm(b);
% Plotting the original matrix result
semilogy(res, 'm-.','LineWidth',2);
hold on

[A, b] = data_prep(1, 1); %=======================================================
[x, res] = custom_conjgrad(A, b, b);
result = norm(A*x -b)/norm(b);
% Plotting the matrix + squares results
semilogy(res, 'b--','LineWidth',2);
hold on

[A, b] = data_prep(1, 2); %=======================================================
[x, res] = custom_conjgrad(A, b, b);
result = norm(A*x -b)/norm(b);
% Plotting the matrix + colwise multiplication results
semilogy(res, 'g:','LineWidth',2);
hold on

[A, b] = data_prep(1, 0); %=======================================================
[x, res] = custom_conjgrad(A, b, b);
result = norm(A*x -b)/norm(b);
% Plotting the matrix + colwise multiplication results
semilogy(res, 'r-','LineWidth',2);
hold on

title('Conjugate Gradient')
xlabel('Dimension of X')

ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$X ^{(1491 \times 9)}$', '$\hat{X}^{(1491\times 18)}$', '$\hat{X} ^{(1491\times 45)}$', '$\hat{X} ^{(1491\times 54)}$','fontsize',16, 'Interpreter', 'Latex')