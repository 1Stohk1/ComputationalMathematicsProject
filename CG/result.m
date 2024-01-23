addpath(genpath([fileparts(pwd), filesep]));

compute_result(-1)
compute_result(1)
compute_result(2)
compute_result(0)


function compute_result(dim)
    [sym_A, sym_b] = data_prep(1, dim); 
    [A, b] = data_prep(0, dim);
    x_star = sym_A\sym_b;

    time = 0;
    iter = 1000;
    for i=1:iter
        [x, ~, ~, t] = custom_conjgrad(sym_A, sym_b, sym_b, 1e-16);
        time = time + t(end);
    end
    time = time / iter;

    result = norm(A*x -b)/norm(b);
    diff_star = norm(x-x_star)/norm(x_star);
    nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
    
    fprintf('A dimension == %.f\n', size(A,2));
    fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
    fprintf('Latex Table \n%d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, time);
end