addpath(genpath([fileparts(pwd), filesep]));

% compute_plot(-1, 'm-.','#ad00ad');
% compute_plot(1, 'b--', '#0000bf');
% compute_plot(2, 'g:', '#00c900');
compute_plot(0, 'r-', '#b50000');

xlim([0 54])
ylim([1e-12 1e3])
title('Conjugate Gradient')
xlabel('Iterations')
ylabel('\boldmath$ ||e_{(i)}||$', 'Interpreter', 'Latex')

% ylabel('\boldmath$ ||e_{(i)}||_A \leq \left( \frac{\sqrt{\kappa} - 1} {\sqrt{\kappa} + 1} \right)^i ||e_{(0)}||_A$', 'Interpreter', 'Latex')
legend('$X ^{(1491 \times 9)}$', '', ...
             '$\hat{X}^{(1491\times 18)}$', '',  ...
             '$\hat{X} ^{(1491\times 45)}$', '',  ...
             '$\hat{X} ^{(1491\times 54)}$', '',  ...
             'fontsize',16,  ...
             'Interpreter', 'Latex',  ...
             'Location', 'southeast')

function compute_plot(dim, str1, str2)
[A, b] = data_prep(1, dim); 
[~, res, ~, ~, c] = custom_conjgrad(A, b, zeros(size(b)), 1e-14, size(A,2), 1);
semilogy(res, str1,'LineWidth', 2);
hold on 
semilogy(c, "Color", str2, 'LineWidth', 2);
end