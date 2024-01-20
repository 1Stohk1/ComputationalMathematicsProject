addpath(genpath([fileparts(pwd), filesep]));

rng(42)
square = 100;
while true
  A = rand(square);
  if rank(A) == square; break; end    %will be true nearly all the time
end
b = rand(square, 1);

b = A'*b;
A = A'*A;
times = zeros(1, square);
for j= 1:1000
    [x, res, iter, ex_time] = custom_conjgrad(A, b, zeros(size(b)), 0);
    times = times + ex_time;
end
% Plotting the matrix + colwise multiplication results
times = times*1000
plot(times, 'r-','LineWidth',2);
xlim([0 square])
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