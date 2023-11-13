addpath(genpath([fileparts(pwd), filesep]));
format long

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

disp('ORIGINAL MATRIX')
va=0;
% Computing the average time for each run
for i=1:1000
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    va = va + toc;
end
x_star =  A\b;
norm(x-x_star)
norm(2*A'*A*x-2*A'*b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1000

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
for i=1:1000
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    va = va + toc;
end
x_star =  A\b;
norm(x-x_star)
norm(2*A'*A*x-2*A'*b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1000

disp('SQUARED +COMBINATION MATRIX')
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
va=0;
% Computing the average time for each run
for i=1:1000
    tic
    [Q, R, x] = custom_opt_HQR(A, b);
    va = va + toc;
end
x_star =  A\b;
norm(x-x_star)
norm(2*A'*A*x-2*A'*b)
norm(A*x -b)/norm(b)
disp('The average time in 1000 runs is:')
va/1000