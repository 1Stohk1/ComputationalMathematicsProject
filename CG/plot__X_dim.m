addpath(genpath([fileparts(pwd), filesep]));

compute_plot(-1, 'm-.');
compute_plot(1, 'b--');
compute_plot(2, 'g:');
compute_plot(0, 'r-');

xlim([0, 54]);
title('Conjugate Gradient')
xlabel('Iterations')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$X ^{(1491 \times 9)}$', ...
             '$\hat{X}^{(1491\times 18)}$',  ...
             '$\hat{X} ^{(1491\times 45)}$',  ...
             '$\hat{X} ^{(1491\times 54)}$', ...
             'fontsize',16, 'Interpreter', 'Latex')
    
function compute_plot(mod, str1)
[A, b] = data_prep(1, mod);
[or_A, or_b] = data_prep(0, mod);
[~, res] = custom_conjgrad(A, b, b, 1e-14, size(A, 2), 0, or_A, or_b);
semilogy(res, str1,'LineWidth', 2);
hold on 
end