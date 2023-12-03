addpath(genpath([fileparts(pwd), filesep]));


[A, b] = data_prep(1, 2); %=======================================================
iter = 1000;
precision  = 1e-8;
time = 0;
for i=1:1000
    tic
    [x, res, atIter] = custom_conjgrad(A, b, b, precision, iter);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f = norm(2*A*x-2*b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

precision  = 1e-10;
time = 0;
for i=1:1000
    tic
    [x, res, atIter] = custom_conjgrad(A, b, b, precision, iter);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f = norm(2*A*x-2*b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);


precision  = 1e-12;
time = 0;
for i=1:1000
    tic
    [x, res, atIter] = custom_conjgrad(A, b, b, precision, iter);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f = norm(2*A*x-2*b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

precision  = 1e-14;
time = 0;
for i=1:1000
    tic
    [x, res, atIter] = custom_conjgrad(A, b, b, precision, iter);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f = norm(2*A*x-2*b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);


precision  = 1e-16;
time = 0;
for i=1:1000
    tic
    [x, res, atIter] = custom_conjgrad(A, b, b, precision, 84);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f = norm(2*A*x-2*b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);