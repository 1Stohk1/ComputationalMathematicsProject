addpath(genpath([fileparts(pwd), filesep]));
format long
fprintf('\n')

[A, b] = data_prep(0, -1); %=======================================================
[Q, R] = custom_HQR(A);
x = R\Q*b;
x_star = A\b;
diff_star = norm(x-x_star)/norm(x_star);
result = norm(A*x-b)/norm(b);
other_res = norm(R*x-Q*b)/norm(b);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('A Size: %d, Result: %d, Other Result: %d, Upper: %d, Lower: %d, Diff_star: %d\n', size(A,2), result, other_res, upper, lower, diff_star)
fprintf(' %d & %d \n', upper, lower)
[Q, R, x] = custom_opt_HQR(A, b);
result = norm(A*x-b)/norm(b);
diff_star = norm(x-x_star);
fprintf('Optimized QR, Result: %d, Diff_star: %d\n', result, diff_star)

[A, b] = data_prep(0, 1); %=======================================================
[Q, R] = custom_HQR(A);
x = R\Q*b;
x_star = A\b;
diff_star = norm(x-x_star);
result = norm(A*x-b)/norm(b);
other_res = norm(R*x-Q*b)/norm(b);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('A Size: %d, Result: %d, Other Result: %d, Upper: %d, Lower: %d, Diff_star: %d\n', size(A,2), result, other_res, upper, lower, diff_star)
fprintf(' %d & %d \n', upper, lower)
[Q, R, x] = custom_opt_HQR(A, b);
result = norm(A*x-b)/norm(b);
diff_star = norm(x-x_star);
fprintf('Optimized QR, Result: %d, Diff_star: %d\n', result, diff_star)

[A, b] = data_prep(0, 2); %=======================================================
[Q, R] = custom_HQR(A);
x = R\Q*b;
x_star = A\b;
diff_star = norm(x-x_star);
result = norm(A*x-b)/norm(b);
other_res = norm(R*x-Q*b)/norm(b);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('A Size: %d, Result: %d, Other Result: %d, Upper: %d, Lower: %d, Diff_star: %d\n', size(A,2), result, other_res, upper, lower, diff_star)
fprintf(' %d & %d \n', upper, lower)
[Q, R, x] = custom_opt_HQR(A, b);
result = norm(A*x-b)/norm(b);
diff_star = norm(x-x_star);
fprintf('Optimized QR, Result: %d, Diff_star: %d\n', result, diff_star)

[A, b] = data_prep(0, 0); %=======================================================
[Q, R] = custom_HQR(A);
x = R\Q*b;
x_star = A\b;
diff_star = norm(x-x_star);
result = norm(A*x-b)/norm(b);
other_res = norm(R*x-Q*b)/norm(b);
Rx = R*x;
Qb = Q*b;
upper = norm(Rx(1:size(A,2),:)-Qb(1:size(A, 2),:))/norm(b);
lower = norm(Qb(size(A, 2)+1:size(Q,1),:))/norm(b);
fprintf('A Size: %d, Result: %d, Other Result: %d, Upper: %d, Lower: %d, Diff_star: %d\n', size(A,2), result, other_res, upper, lower, diff_star)
fprintf(' %d & %d \n', upper, lower)
[Q, R, x] = custom_opt_HQR(A, b);
result = norm(A*x-b)/norm(b);
diff_star = norm(x-x_star);
fprintf('Optimized QR, Result: %d, Diff_star: %d\n', result, diff_star)