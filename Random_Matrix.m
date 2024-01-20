rng(42)

matrixColumn = 200;
matrixRow = 2000;
while true
  A = rand(matrixRow, matrixColumn);
  if rank(A) == matrixColumn; break; end    %will be true nearly all the time
end
b = rand(matrixRow, 4);
x_star = A\b;
result = norm(A*x_star-b)/norm(b);

time_qr = 0;
for i=1:100
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/100;
QResult = norm(A*x-b)/norm(b);
diff_qr = norm(x_star-x)/norm(x_star);
nabla_qr = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

symm_b = A'*b;
symm_A = A'*A;
time_cg = 0;
for i=1:100
    [x, res, iter, ex_time] = custom_conjgrad(symm_A, symm_b, symm_b, 5e-15);
    time_cg = time_cg + ex_time(iter);
end
time_cg = time_cg/100;
cg_result = norm(A*x-b)/norm(b);
diff_cg = norm(x_star-x)/norm(x_star);
nabla_cg = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

fprintf('Lib Result: %d \n  QR Result: %d \n CG Result: %d\n', result, QResult, cg_result)
fprintf('%d & %d  & %d  & %d \\ %d & %d  & %d & %d \n', diff_qr, nabla_qr, QResult, time_qr, diff_cg, nabla_cg, cg_result, time_cg)


square = 800;
while true
  A = rand(square);
  if rank(A) == square; break; end    %will be true nearly all the time
end
b = rand(square, 2);
x_star = A\b;
result = norm(A*x_star-b)/norm(b);

time_qr = 0;
for i=1:10
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/10;
QResult = norm(A*x-b)/norm(b);
diff_qr = norm(x_star-x)/norm(x_star);
nabla_qr = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

symm_b = A'*b;
symm_A = A'*A;
x_star = A\b;
time_cg = 0;
for i=1:10
    [x, res, iter, ex_time] = custom_conjgrad(symm_A, symm_b, zeros(size(symm_b)), 5e-15, 1500);
    time_cg = time_cg + ex_time(iter);
end
time_cg = time_cg/10;
cg_result = norm(A*x-b)/norm(b);
diff_cg = norm(x_star-x)/norm(x_star);
nabla_cg = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);


fprintf('Lib Result: %d \n  QR Result: %d \n CG Result: %d\n', result, QResult, cg_result)
fprintf('%d & %d  & %d  & %d \\ %d & %d  & %d & %d \n', diff_qr, nabla_qr, QResult, time_qr, diff_cg, nabla_cg, cg_result, time_cg)
