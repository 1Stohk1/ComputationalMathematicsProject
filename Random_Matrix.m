matrixColumn = 800;
matrixRow = 2000;
rng(42)
while true
  A = rand(matrixRow, matrixColumn);
  if rank(A) == matrixColumn; break; end    %will be true nearly all the time
end

b = rand(matrixRow, 4);
x_star = A\b;
result = norm(A*x_star-b)/norm(b);

[Q, R, x] = custom_opt_HQR(A, b);
QResult = norm(A*x-b)/norm(b);
diff_qr = norm(x_star-x)/norm(x_star);
nabla_qr = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

b = A'*b;
A = A'*A;
% x_star = A\b;
[x, res, iter, ex_time] = custom_conjgrad(A, b, b, 1e-14);
cg_result = norm(A*x-b)/norm(b);
diff_cg = norm(x_star-x)/norm(x_star);
nabla_cg = norm(2*(A)*x-2*b)/norm(A);


fprintf('Lib Result: %d \n  QR Result: %d \n CG Result: %d\n', result, QResult, cg_result)
fprintf('%d & %d  & %d \\ %d & %d  & %d \n', diff_qr, nabla_qr, QResult, diff_cg, nabla_cg, cg_result)

square = 1400;
while true
  A = rand(square);
  if rank(A) == square; break; end    %will be true nearly all the time
end

b = rand(square, 4);
x_star = A\b;
result = norm(A*x_star-b)/norm(b);

[Q, R, x] = custom_opt_HQR(A, b);
QResult = norm(A*x-b)/norm(b);
diff_qr = norm(x_star-x)/norm(x_star);
nabla_qr = norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);

b = A'*b;
A = A'*A;
x_star = A\b;
[x, res, iter, ex_time] = custom_conjgrad(A, b, b, 1e-14);
cg_result = norm(A*x-b)/norm(b);
diff_cg = norm(x_star-x)/norm(x_star);
nabla_cg = norm(2*(A)*x-2*b)/norm(A);


fprintf('Lib Result: %d \n  QR Result: %d \n CG Result: %d\n', result, QResult, cg_result)
fprintf('%d & %d  & %d \\ %d & %d  & %d \n', diff_qr, nabla_qr, QResult, diff_cg, nabla_cg, cg_result)