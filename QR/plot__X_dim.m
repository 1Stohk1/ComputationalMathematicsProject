addpath(genpath([fileparts(pwd), filesep]));

compute_plot(2, 'g-');
compute_plot(0, 'r-');

xlim([0 54])
    
title('Householder QR')
xlabel('Dimension of X')
ylabel('\boldmath$ ||\hat{X} w-y||\slash||y||$', 'Interpreter', 'Latex')
legend('$\hat{X} ^{(1491\times 45)}$', '', '$\hat{X} ^{(1491\times 54)}$', '$||Q^T_2y||\slash||y||$','fontsize',16, 'Interpreter', 'Latex')

function compute_plot(dim, line)
    [A, b] = data_prep(0, dim);
    [m, n] = size(A);
    
    atIter = 0;
    results = zeros(1, n); 
    residuals = zeros(1, n);
    for k = 1:n
        atIter = atIter + 1;
        [Q1, R1] = custom_HQR(A(:, 1:k));
        x = R1\(Q1*b);
        results(atIter) = norm(R1*x-Q1*b)/norm(b);
        Qb = Q1*b;
        residuals(atIter) = norm(Qb(k+1:m,:))/norm(b);
    end

    semilogy(results, line,'LineWidth',2);
    hold on
    semilogy(residuals, 'b.', 'LineWidth', 1.5);
end



