addpath(genpath([fileparts(pwd), filesep]));

[time_qr, result_qr, diff_star_qr, nabla_f_qr] = computing_qr(-1);
[time_cg, result_cg, diff_star_cg, nabla_f_cg] = computing_cg(-1, diff_star_qr);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[time_qr, result_qr, diff_star_qr, nabla_f_qr] = computing_qr(1);
[time_cg, result_cg, diff_star_cg, nabla_f_cg] = computing_cg(1, diff_star_qr);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[time_qr, result_qr, diff_star_qr, nabla_f_qr] = computing_qr(2);
[time_cg, result_cg, diff_star_cg, nabla_f_cg] = computing_cg(2, diff_star_qr);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)

[time_qr, result_qr, diff_star_qr, nabla_f_qr] = computing_qr(0);
[time_cg, result_cg, diff_star_cg, nabla_f_cg] = computing_cg(0, diff_star_qr);
fprintf('makecell{%d \\ %d} & makecell{%d \\ %d}  & makecell{%d \\ %d}  & makecell{%d \\ %d}  \n', diff_star_qr, diff_star_cg,  nabla_f_qr, nabla_f_cg, result_qr, result_cg, time_qr, time_cg)


function[time, result, diff_star, nabla] = computing_qr(dim, print)
if nargin < 2
    print = 0;
end
    [A, b] = data_prep(0, dim);
    x_star =  A\b;

    time = 0;
    for i=1:1000
        tic
        [~, ~, x] = custom_opt_HQR(A, b);
        time = time + toc;
    end
    time = time/1000;
    
    result = norm(A*x -b)/norm(b);
    diff_star = norm(x-x_star)/norm(x_star);
    nabla =norm(2*(A'*A)*x-2*A'*b)/norm(A'*A);
    
if print == 1
    fprintf('Result: %d; \nDifference: %d; \nNabla: %d', result, diff_star, nabla);
end
end

function[time, result, diff_star, nabla] = computing_cg(dim, threshold, print)
if nargin < 3
    print = 0;
end
[sym_A, sym_b] = data_prep(1, dim);
[A, b] = data_prep(0, dim);
x_star =  A\b;

time = 0;
for i=1:1000
    [x_cg, ~, it, t] = custom_conjgrad(sym_A, sym_b, sym_b, threshold*5, size(A, 2)*4, 0);
    time = time + t(it);
end
time = time/1000;

result = norm(A*x_cg -b)/norm(b);
diff_star= norm(x_cg-x_star)/norm(x_star);
nabla= norm(2*(A'*A)*x_cg-2*A'*b)/norm(A'*A);
if print == 1
    fprintf('Result: %d; \nDifference: %d; \nNabla: %d', result, diff_star, nabla);
end
end
