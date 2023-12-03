addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(0, 0); %=======================================================
tic
[m, n] = size(A);
atIter = 0;
residuals = zeros(1, n); 
for k = 1:n
    atIter = atIter + 1;
    [Q, R, x] = custom_opt_HQR(A(:, 1:k), b);
    residuals(atIter) = norm(A(:, 1:k)*x -b)/norm(b);
end
toc
plot(residuals, 'LineWidth', 2);

title('Householder QR')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
% legend('Original A', 'A & Squares Cols', 'A & Colwise Mul', 'A & Squares & Colwise')
