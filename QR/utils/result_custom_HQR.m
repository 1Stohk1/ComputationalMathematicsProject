addpath(genpath([fileparts(pwd), filesep]));
format long
[A, b] = data_prep(0,0);
[Q,R] = custom_HQR(A);

x_qr = R\(Q*b);
x_star =  A\b;
result = norm(A*x_qr-b)/norm(b);
result_star = norm(A*x_star-b)/norm(b);
fprintf('Result Star %d, Result QR %d \n',result_star, result)

[A, b] = data_prep(1, 0); %=======================================================
[x, r, i, t] = custom_conjgrad(A, b, b, 1e-16);
[A, b] = data_prep(0,0);
result = norm(A*x-b)/norm(b);
fprintf('Result Star %d, Result CG %d \n',result_star, result)