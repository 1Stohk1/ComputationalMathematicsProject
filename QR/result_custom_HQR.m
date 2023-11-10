addpath(genpath([fileparts(pwd), filesep]));
format long

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

disp('ORIGINAL MATRIX')
va=0;
% Computing the average time for each run
for i=1:1
    tic
    [Q, R] = custom_HQR(A);
    va = toc;
end
x = R\(Q\b);
x_star =  A\b;
norm(x-x_star)
norm(A*x -b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1

disp('SQUARED MATRIX')
v = 1:size(A,2);
c = nchoosek(v,2);
for i = 1:size(A,2)
    NewCol = abs(A(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
va=0;
% Computing the average time for each run
for i=1:1
    tic
    [Q, R] = custom_HQR(A);
    va = toc;
end
x = R\(Q\b);
x_star =  A\b;
norm(x-x_star)
norm(A*x -b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1

disp('SQUARED +COMBINATION MATRIX')
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
clear NewCol i temp v c
va=0;
% Computing the average time for each run
for i=1:1
    tic
    [Q, R] = custom_HQR(A);
    va = toc;
end
x = R\(Q\b);
x_star =  A\b;
norm(x-x_star)
norm(A*x -b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1

% tic
% [Q, R] = custom_HQR(A);
% x = R\(Q\b);
% toc
% norm(A*x -b)/norm(b)
% x_star =  A\b;
% norm(x-x_star)

% 
% tic
% [Q, R] = custom_HQR(A);
% x = R\(Q\b);
% toc
% norm(A*x -b)/norm(b)
% 
% tic
% [Q, R] = qr(A);
% x = R\(Q\b);
% toc
% norm(A*x -b)/norm(b)
% 
% tic
% x = A\b;
% toc
% norm(A*x - b)