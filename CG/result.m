addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(1, -1); %=======================================================
time = 0;
for i=1:1000
    tic
    [x] = custom_conjgrad(A, b, b, 1e-16);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*A*x-2*b)/norm(A);
fprintf('A dimension == %.f\n', size(A,1));
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d \n\n', size(A,1), diff_star, nabla_f, result, time);

[A, b] = data_prep(1, 1); %=======================================================
time = 0;
for i=1:1000
    tic
    [x] = custom_conjgrad(A, b, b, 1e-16);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*A*x-2*b)/norm(A);
fprintf('A dimension == %.f\n', size(A,1));
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d \n\n', size(A,1), diff_star, nabla_f, result, time);

[A, b] = data_prep(1, 2); %=======================================================
time = 0;
for i=1:1000
    tic
    [x] = custom_conjgrad(A, b, b, 1e-16);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*A*x-2*b)/norm(A);
fprintf('A dimension == %.f\n', size(A,1));
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d \n\n', size(A,1), diff_star, nabla_f, result, time);

[A, b] = data_prep(1, 0); %=======================================================
time = 0;
for i=1:1000
    tic
    [x] = custom_conjgrad(A, b, b, 1e-16);
    time = time + toc;
end
x_star =  A\b;
time = time/1000;
result = norm(A*x -b)/norm(b);
diff_star = norm(x-x_star)/norm(x_star);
nabla_f = norm(2*A*x-2*b)/norm(A);
fprintf('A dimension == %.f\n', size(A,1));
fprintf('Difference: %d\nNabla: %d\nResult: %d\nTime: %d\n', diff_star, nabla_f, result, time);
fprintf('%d & %d & %d & %d & %d \n\n', size(A,1), diff_star, nabla_f, result, time);


% LIBRARY =================================================================
% [ma,na]=size(A);
% [mb,nb]=size(b);
% afun=@(x)  reshape(A*reshape(x,na,[]),[],1);
% x=pcg(afun,b(:), 1e-8);
% x=reshape(x,na,nb);
%
% result = norm(A*x -b)/norm(b);
% diff_star = norm(x-x_star);
% nabla_f = norm(2*A*x-2*b);
% fprintf('LIBRARY \nDIFFERENCE: %d\nNABLA: %d\nRESULT: %d\n', diff_star, nabla_f, result);


