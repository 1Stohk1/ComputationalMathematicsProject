addpath(genpath([fileparts(pwd), filesep]));


[A, b] = data_prep(1, 2); %=======================================================
[x, res, iter, ex_time] = custom_conjgrad(A, b, zeros(size(b)), 0, 100);
result = norm(A*x -b)/norm(b);
% Plotting the matrix + colwise multiplication results
plot(ex_time, 'r-','LineWidth',2);
hold on

% [ma,na]=size(A);
% [mb,nb]=size(b);
% afun=@(x)  reshape(A*reshape(x,na,[]),[],1);
% [X,FLAG,RELRES,ITER,RESVEC]=pcg(afun,b(:), 1e-16, 1000);
% x=reshape(x,na,nb);
% % Plotting the matrix + colwise multiplication results
% RESVEC = RESVEC/norm(b);
% semilogy(RESVEC, 'r-.','LineWidth',2);
% hold on

title('Conjugate Gradient (Time)')
xlabel('Iterations')
ylabel('Ex. Time ($\mu$)', 'Interpreter', 'Latex')
% legend('$\hat{X} ^{(1491\times 54)}$', 'Interpreter', 'Latex')