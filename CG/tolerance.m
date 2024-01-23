addpath(genpath([fileparts(pwd), filesep]));

tolerance_computation(1e-2);
tolerance_computation(1e-3);
tolerance_computation(1e-4);
tolerance_computation(1e-14);
tolerance_computation(1e-16);

function tolerance_computation(precision)
    [sym_A, sym_b] = data_prep(1, 0);
    [A, b] = data_prep(0, 0); 
    x_star =  A\b;
    
    time = 0;
    iter = 1000;
    for i=1:iter
        [x, ~, atIter, t] = custom_conjgrad(sym_A, sym_b, sym_b, precision, 126);
        time = time + t(atIter);
    end
    time = time/iter;

    result = norm(A*x -b)/norm(b);
    diff_star = norm(x-x_star)/norm(x_star);
    nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

    fprintf('A dimension == %.f\n', size(A,2));
    fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
    fprintf('Latex Table \n%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);
end