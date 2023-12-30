addpath(genpath([fileparts(pwd), filesep]));
format long

[A, b] = data_prep(0, -1); %=======================================================
time = 0;
for i=1:100
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time = time + toc;
end
time = time/100;
result = norm(A*x -b)/norm(b);
x_star =  A\b;
diff_star = norm(x-x_star)/norm(x_star);
result_star = norm(A*x_star -b)/norm(b);
nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('A dimension == %.f\n', size(A,2));

[Q, R] = custom_HQR(A);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nUpper Block: %d\nTime: %d\n', diff_star, nabla_f, result, upper, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, upper, time);


[A, b] = data_prep(0, 1); %=======================================================
time = 0;
for i=1:100
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time = time + toc;
end
time = time/100;
result = norm(A*x -b)/norm(b);
x_star =  A\b;
diff_star = norm(x-x_star)/norm(x_star);
result_star = norm(A*x_star -b)/norm(b);
nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('A dimension == %.f\n', size(A,2));

[Q, R] = custom_HQR(A);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nUpper Block: %d\nTime: %d\n', diff_star, nabla_f, result, upper, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, upper, time);


[A, b] = data_prep(0, 2); %=======================================================
time = 0;
for i=1:100
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time = time + toc;
end
time = time/100;
result = norm(A*x -b)/norm(b);
x_star =  A\b;
diff_star = norm(x-x_star)/norm(x_star);
result_star = norm(A*x_star -b)/norm(b);
nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('A dimension == %.f\n', size(A,2));

[Q, R] = custom_HQR(A);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nUpper Block: %d\nTime: %d\n', diff_star, nabla_f, result, upper, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, upper, time);


[A, b] = data_prep(0, 0); %=======================================================
time = 0;
for i=1:100
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time = time + toc;
end
time = time/100;
result = norm(A*x -b)/norm(b);
x_star =  A\b;
diff_star = norm(x-x_star)/norm(x_star);
result_star = norm(A*x_star -b)/norm(b);
nabla_f =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
fprintf('A dimension == %.f\n', size(A,2));

[Q, R] = custom_HQR(A);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('Difference: %d\nNabla: %d\nResult: %d\nUpper Block: %d\nTime: %d\n', diff_star, nabla_f, result, upper, time);
fprintf('%d & %d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, upper, time);