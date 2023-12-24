addpath(genpath([fileparts(pwd), filesep]));
format long

[A, b] = data_prep(0, -1); %=======================================================
[Q, R, x] = custom_opt_HQR(A, b);
x_star =  A\b;
result = norm(A*x-b);
diff_star = norm(x-x_star);
result_star = norm(A*x_star-b);
fprintf('A dimension == %.f\n', size(A,2));
fprintf('Diff Star: %d\nDiff b: %d\nDiff b_star: %d\n', diff_star, result, result_star);


[A, b] = data_prep(0, 0); %=======================================================
[Q, R, x] = custom_opt_HQR(A, b);
x_star =  A\b;
result = norm(A*x-b);
diff_star = norm(x-x_star);
result_star = norm(A*x_star-b);
fprintf('A dimension == %.f\n', size(A,2));
fprintf('Diff Star: %d\nDiff b: %d\nDiff b_star: %d\n', diff_star, result, result_star);
fprintf('%d & %d  \n\n', size(A,2), result);