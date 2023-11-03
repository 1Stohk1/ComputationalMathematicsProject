addpath(genpath([fileparts(pwd), filesep]));

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);
symm_b = A'*b;
symm_A = A'*A;
[x, res] = custom_conjgrad(symm_A, symm_b);
result = symm_A*x -symm_b;
norm(result)
% Plotting the original matrix result
plot(res, 'r-','LineWidth',2);
hold on
clear x res symm_A symm_b

A1 = A;
v = 1:size(A1,2);
c = nchoosek(v,2);
for i = 1:size(A1,2)
    NewCol = abs(A1(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A1 = [A1 NewCol];
end
symm_b = A1'*b;
symm_A = A1'*A1;
[x, res] = custom_conjgrad(symm_A, symm_b);
result = symm_A*x -symm_b;
norm(result)
% Plotting the matrix + squares results
plot(res, 'b--','LineWidth',2);
hold on
clear x res symm_A symm_b A1

A2 = A;
v = 1:size(A2,2);
c = nchoosek(v,2);
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A2(:,j).*A2(:,k));
    NewCol = NewCol/norm(NewCol);
    A2 = [A2 NewCol];
end
clear NewCol i temp v c
symm_b = A2'*b;
symm_A = A2'*A2;
[x, res] = custom_conjgrad(symm_A, symm_b);
result = symm_A*x -symm_b;
norm(result)
% Plotting the matrix + colwise multiplication results
plot(res, 'g:','LineWidth',2);
hold on
clear x res symm_A symm_b A2

A3 = A;
v = 1:size(A3,2);
c = nchoosek(v,2);
for i = 1:size(A3,2)
    NewCol = abs(A3(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A3 = [A3 NewCol];
end
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A3(:,j).*A3(:,k));
    NewCol = NewCol/norm(NewCol);
    A3 = [A3 NewCol];
end
clear NewCol i temp v c
symm_b = A3'*b;
symm_A = A3'*A3;
[x, res] = custom_conjgrad(symm_A, symm_b);
result = symm_A*x -symm_b;
norm(result)
% Plotting the matrix + colwise multiplication results
plot(res, 'm-.','LineWidth',2);
hold on
clear x res symm_A symm_b A3

title('Conjugate Gradient')
xlabel('Iteration')
ylabel('Loss')
legend('Original A', 'A & Squares Cols', 'A & Colwise Mul', 'A & Squares & Colwise')