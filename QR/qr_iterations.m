addpath(genpath([fileparts(pwd), filesep]));

rng(42)
square = 100;
while true
  A = rand(square);
  if rank(A) == square; break; end    %will be true nearly all the time
end
b = rand(square, 1);

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
% plot(results, 'LineWidth', 2);
% hold on
% plot(residuals, 'r--', 'LineWidth', 1.5);
% 
% title('Householder QR')
% xlabel('Dimension of X')
% ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
% legend('$||Rw-Q^Ty||\slash||y||$', '$||Q^T_2y||\slash||y||$', 'fontsize',16, 'Interpreter', 'Latex')
% 
% clf;

times = times  * 100;

plot(times, 'r-','LineWidth',2);
xlim([0 square])
title('QR Factorization (Time)')
xlabel('Iterations')
ylabel('Ex. Time ($\mu$)', 'Interpreter', 'Latex')