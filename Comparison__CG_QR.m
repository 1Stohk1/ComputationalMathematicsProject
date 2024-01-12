addpath(genpath([fileparts(pwd), filesep]));

[A, b] = data_prep(0, -1); %=======================================================
fprintf('A dimension == %.f\n', size(A,2));
x_star =  A\b;
result_star = norm(A*x_star -b)/norm(b);
time_qr = 0;
for i=1:1000
    tic
    [Q, R, x_qr] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/1000;
result_qr = norm(A*x_qr -b)/norm(b);
diff_star_qr = norm(x_qr-x_star)/norm(x_star);
nabla_f_qr =norm(2*(A'*A)*x_qr-2*A'*b)/norm(A'*A);
[A, b] = data_prep(1, -1); %=======================
x_star =  A\b;
time_cg = 0;
for i=1:1000
    [x_cg, r, i, t] = custom_conjgrad(A, b, b, diff_star_qr*2, 300);
    time_cg = time_cg + t(i);
end
time_cg = time_cg/1000;
result_cg = norm(A*x_cg -b)/norm(b);
diff_star_cg = norm(x_cg-x_star)/norm(x_star);
nabla_f_cg = norm(2*A*x_cg-2*b)/norm(A);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[A, b] = data_prep(0, 1); %=======================================================
fprintf('A dimension == %.f\n', size(A,2));
x_star =  A\b;
result_star = norm(A*x_star -b)/norm(b);
time_qr = 0;
for i=1:1000
    tic
    [Q, R, x_qr] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/1000;
result_qr = norm(A*x_qr -b)/norm(b);
diff_star_qr = norm(x_qr-x_star)/norm(x_star);
nabla_f_qr =norm(2*(A'*A)*x_qr-2*A'*b)/norm(A'*A);
[A, b] = data_prep(1, 1); %=======================
x_star =  A\b;
time_cg = 0;
for i=1:1000
    [x_cg, r, i, t] = custom_conjgrad(A, b, b, diff_star_qr*2, 300);
    time_cg = time_cg + t(i);
end
time_cg = time_cg/1000;
result_cg = norm(A*x_cg -b)/norm(b);
diff_star_cg = norm(x_cg-x_star)/norm(x_star);
nabla_f_cg = norm(2*A*x_cg-2*b)/norm(A);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[A, b] = data_prep(0, 2); %=======================================================
fprintf('A dimension == %.f\n', size(A,2));
x_star =  A\b;
result_star = norm(A*x_star -b)/norm(b);
time_qr = 0;
for i=1:1000
    tic
    [Q, R, x_qr] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/1000;
result_qr = norm(A*x_qr -b)/norm(b);
diff_star_qr = norm(x_qr-x_star)/norm(x_star);
nabla_f_qr =norm(2*(A'*A)*x_qr-2*A'*b)/norm(A'*A);
[A, b] = data_prep(1, 2); %=======================
x_star =  A\b;
time_cg = 0;
for i=1:1000
    [x_cg, r, i, t] = custom_conjgrad(A, b, b, diff_star_qr*2, 300);
    time_cg = time_cg + t(i);
end
time_cg = time_cg/1000;
result_cg = norm(A*x_cg -b)/norm(b);
diff_star_cg = norm(x_cg-x_star)/norm(x_star);
nabla_f_cg = norm(2*A*x_cg-2*b)/norm(A);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[A, b] = data_prep(0, 0); %=======================================================
fprintf('A dimension == %.f\n', size(A,2));
x_star =  A\b;
result_star = norm(A*x_star -b)/norm(b);
time_qr = 0;
for i=1:1000
    tic
    [Q, R, x_qr] = custom_opt_HQR(A, b);
    time_qr = time_qr + toc;
end
time_qr = time_qr/1000;
result_qr = norm(A*x_qr -b)/norm(b);
diff_star_qr = norm(x_qr-x_star)/norm(x_star);
nabla_f_qr =norm(2*(A'*A)*x_qr-2*A'*b)/norm(A'*A);
[A, b] = data_prep(1, 0); %=======================
x_star =  A\b;
time_cg = 0;
for i=1:1000
    [x_cg, r, i, t] = custom_conjgrad(A, b, b, diff_star_qr*2, 300);
    time_cg = time_cg + t(i);
end
time_cg = time_cg/1000;
result_cg = norm(A*x_cg -b)/norm(b);
diff_star_cg = norm(x_cg-x_star)/norm(x_star);
nabla_f_cg = norm(2*A*x_cg-2*b)/norm(A);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)