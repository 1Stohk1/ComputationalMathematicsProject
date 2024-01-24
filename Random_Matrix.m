[A, b] = system_generator(200, 2000, 4);

[time_qr, result_qr, diff_qr, nabla_qr] = computing_qr(A, b, 10); % It is impossible to average more than this
[time_cg, result_cg, diff_cg, nabla_cg] = computing_cg(A, b, 10); % It is impossible to average more than this

fprintf('Latex Table\n');
fprintf('%d & %d  & %d  & %d \\ %d & %d  & %d & %d \n', ...
            diff_qr, nabla_qr, result_qr, time_qr, ...
            diff_cg, nabla_cg, result_cg, time_cg)

[A, b] = system_generator(800, 800, 4);

[time_qr, result_qr, diff_qr, nabla_qr] = computing_qr(A, b, 5); % It is impossible to average more than this
[time_cg, result_cg, diff_cg, nabla_cg] = computing_cg(A, b, 5); % It is impossible to average more than this

fprintf('Latex Table\n');
fprintf('%d & %d  & %d  & %d \\ %d & %d  & %d & %d \n', ...
            diff_qr, nabla_qr, result_qr, time_qr, ...
            diff_cg, nabla_cg, result_cg, time_cg)
        
function[time, result, diff, nabla] = computing_cg(A, b, iters)
    symm_b = A'*b;
    symm_A = A'*A;
    x_star = A\b;
    
    time = 0;
    for i=1:iters
        [x, ~, iter, ex_time] = custom_conjgrad(symm_A, symm_b, symm_b, 5e-15);
        time = time + ex_time(iter);
    end
    time = time/iters;
    result = norm(A*x-b)/norm(b);
    diff = norm(x_star-x)/norm(x_star);
    nabla = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
end

function[time, result, diff, nabla] = computing_qr(A, b, iters)
    x_star = A\b;
    time = 0;
    for i=1:iters
        tic
        [~, ~, x] = custom_opt_HQR(A, b);
        time = time + toc;
    end
    time = time/iters;
    result = norm(A*x-b)/norm(b);
    diff = norm(x_star-x)/norm(x_star);
    nabla = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
end