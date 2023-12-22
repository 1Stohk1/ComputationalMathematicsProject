addpath(genpath([fileparts(pwd), filesep]));
format long

[A, b] = data_prep(0, -1); %=======================================================
time = 0;
for i=1:10
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time = time + toc;
end
x_star =  A\b;
time = time/10;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star);
nabla_f =norm(2*A'*A*x-2*A'*b);
fprintf('A dimension == %.f\n', size(A,2));
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, time);
% 
% [A, b] = data_prep(0, 1); %=======================================================
% time = 0;
% for i=1:1000
%     tic
%     [Q, R, x] = custom_opt_HQR(A, b);
%     time = time + toc;
% end
% x_star =  A\b;
% time = time/1000;
% result = norm(A*x -b)/norm(b);
% diff_star = norm(x-x_star);
% nabla_f =norm(2*A'*A*x-2*A'*b);
% fprintf('A dimension == %.f\n', size(A,2));
% fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
% fprintf('%d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, time);
% 
% [A, b] = data_prep(0, 2); %=======================================================
% time = 0;
% for i=1:1000
%     tic
%     [Q, R, x] = custom_opt_HQR(A, b);
%     time = time + toc;
% end
% x_star =  A\b;
% time = time/1000;
% result = norm(A*x -b)/norm(b);
% diff_star = norm(x-x_star);
% nabla_f =norm(2*A'*A*x-2*A'*b);
% fprintf('A dimension == %.f\n', size(A,2));
% fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
% fprintf('%d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, time);
% 
% [A, b] = data_prep(0, 0); %=======================================================
% time = 0;
% for i=1:1000
%     tic
%     [Q, R, x] = custom_opt_HQR(A, b);
%     time = time + toc;
% end
% x_star =  A\b;
% time = time/1000;
% result = norm(A*x -b)/norm(b);
% diff_star = norm(x-x_star);
% nabla_f =norm(2*A'*A*x-2*A'*b);
% fprintf('A dimension == %.f\n', size(A,2));
% fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
% fprintf('%d & %d & %d & %d & %d \n\n', size(A,2), diff_star, nabla_f, result, time);