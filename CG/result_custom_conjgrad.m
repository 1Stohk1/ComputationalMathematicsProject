addpath(genpath([fileparts(pwd), filesep]));
format long
temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

disp('ORIGINAL MATRIX')
symm_b = A'*b;
symm_A = A'*A;
va=0;
% Computing the average time for each run
for i=1:1000
    tic
    [x, residual] = custom_conjgrad(symm_A, symm_b);
    va = va + toc;
end
[x, residual] = custom_conjgrad(symm_A, symm_b);
x_star =  symm_A\symm_b;
norm(x-x_star)
norm(symm_A*x -symm_b)
norm(symm_A*x -symm_b)/norm(symm_b)
disp('The average time in 1000 runs is:')
va/1000

disp('SQUARED MATRIX')
v = 1:size(A,2);
c = nchoosek(v,2);
for i = 1:v
    NewCol = abs(A(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
symm_b = A'*b;
symm_A = A'*A;
va=0;
% Computing the average time for each run
for i=1:1000
    tic
    [x, residual] = custom_conjgrad(symm_A, symm_b);
    va = va + toc;
end
[x, residual] = custom_conjgrad(symm_A, symm_b);
x_star =  symm_A\symm_b;
norm(x-x_star)
norm(symm_A*x -symm_b)
norm(symm_A*x -symm_b)/norm(symm_b)
disp('The average time in 1000 runs is:')
va/1000

disp('SQUARED + COMBINATION MATRIX')
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
clear NewCol i temp v c
symm_b = A'*b;
symm_A = A'*A;
va=0;
% Computing the average time for each run
for i=1:1000
    tic
    [x, residual] = custom_conjgrad(symm_A, symm_b);
    va = va + toc;
end
[x, residual] = custom_conjgrad(symm_A, symm_b);
x_star =  symm_A\symm_b;
norm(x-x_star)
norm(symm_A*x -symm_b)
norm(symm_A*x -symm_b)/norm(symm_b)
disp('The average time in 1000 runs is:')
va/1000

% tic
% [x, residual] = custom_conjgrad(symm_A, symm_b);
% toc
% norm(symm_A*x -symm_b)
% 
% tic
% [ma,na]=size(symm_A);
% [mb,nb]=size(symm_b);
% afun=@(x)  reshape(symm_A*reshape(x,na,[]),[],1);
% X=pcg(afun,symm_b(:));
% X=reshape(X,na,nb);
% toc
% norm(symm_A*X -symm_b)
% 
% tic
% x = symm_A\symm_b;
% toc
% norm(symm_A*x -symm_b)
