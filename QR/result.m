addpath(genpath([fileparts(pwd), filesep]));

compute_qr(-1);
compute_qr(1);
compute_qr(2);
compute_qr(0);

function compute_qr(dim)
    [A, b] = data_prep(0, dim);
    x_star =  A\b;

    time = 0;
    for i=1:1000
        tic
        [~, ~, x] = custom_opt_HQR(A, b);
        time = time + toc;
    end
    time = time/1000;

    result = norm(A*x -b)/norm(b);
    diff_star = norm(x-x_star)/norm(x_star);
    nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

    [Q, R] = custom_HQR(A);
    Rx = R*x;
    Qb = Q*b;
    upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
    lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);

    fprintf('A dimension == %.f\n', size(A,2));
    fprintf('Difference: %d\nNabla: %d\nResult: %d\nUpper Block: %d\nTime: %d\n', diff_star, nabla_f, result, upper, time);
    fprintf('%d & %d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, upper, time);
end