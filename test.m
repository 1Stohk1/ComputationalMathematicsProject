matrixColumn = 100;
matrixRow = 1000;
rng(42)
while true
  A = rand(matrixRow, matrixColumn);
  if rank(A) == matrixColumn; break; end    %will be true nearly all the time
end

b = rand(matrixRow, 1);
x_star = A\b;
result = norm(A*x_star-b)/norm(b);

[Q, R, x] = custom_opt_HQR(A, b);
ouResult = norm(A*x-b)/norm(b);
diff_qr = norm(x_star-x);

b = A'*b;
A = A'*A;
[x, res, iter, ex_time] = custom_conjgrad(A, b, b, 0, 300);
cg_result = norm(A*x-b)/norm(b);
diff_cg = norm(x_star-x);


fprintf(' %d & %d & %d\n', result, ouResult, cg_result)
fprintf('%d & %d \n', diff_qr, diff_cg)