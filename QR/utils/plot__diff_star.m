addpath(genpath([fileparts(pwd), filesep]));

compute_plot(2, 'g:');
compute_plot(0, 'r-');

xlim([0 54])
title('Householder QR (Difference w)')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$\hat{X} ^{(1491\times 45)}$', '$\hat{X} ^{(1491\times 54)}$' ,'fontsize',16, 'Interpreter', 'Latex')

function compute_plot(dim, line)
    [A, b] = data_prep(0, dim);
    [~, n] = size(A);
    x_star = A\b;
    
    atIter = 0;
    results = zeros(1, n); 
    x_comp = zeros(size(x_star));
    for k = 1:n
        atIter = atIter + 1;
        [~, ~, x] = custom_opt_HQR(A(:, 1:k), b);
        x_comp(1:k,:) = x;
        results(atIter) = norm(x_comp-x_star)/norm(x_star);
    end

    semilogy(results, line,'LineWidth',2);
    hold on
end



