addpath(genpath([fileparts(pwd), filesep]));


[A, b] = data_prep(1, 0); %=======================================================
iter = 1000;
precision  = 1e-2;
time = 0;
for i=1:1000
    [x, res, atIter, t] = custom_conjgrad(A, b, b, precision, iter);
    time = time + t(atIter);
end
[A, b] = data_prep(0, 0); %=======================================================
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

[A, b] = data_prep(1, 0); %=======================================================
precision  = 1e-3;
time = 0;
for i=1:1000
    [x, res, atIter, t] = custom_conjgrad(A, b, b, precision, iter);
    time = time + t(atIter);
end
[A, b] = data_prep(0, 0); %=======================================================
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

[A, b] = data_prep(1, 0); %=======================================================
precision  = 1e-4;
time = 0;
for i=1:1000
    [x, res, atIter, t] = custom_conjgrad(A, b, b, precision, iter);
    time = time + t(atIter);
end
[A, b] = data_prep(0, 0); %=======================================================
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

[A, b] = data_prep(1, 0); %=======================================================
precision  = 1e-8;
time = 0;
for i=1:1000
    [x, res, atIter, t] = custom_conjgrad(A, b, b, precision, iter);
    time = time + t(atIter);
end
[A, b] = data_prep(0, 0); %=======================================================
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);

[A, b] = data_prep(1, 0); %=======================================================
precision  = 1e-12;
time = 0;
for i=1:1000
    [x, res, atIter, t] = custom_conjgrad(A, b, b, precision, 125);
    time = time + t(atIter);
end
[A, b] = data_prep(0, 0); %=======================================================
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', precision, diff_star, nabla_f, result, atIter, time);