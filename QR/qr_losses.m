addpath(genpath([fileparts(pwd), filesep]));

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);
[Q, R] = custom_HQR(A);
fact_A = Q* R;
x = R\(Q\b);
result = fact_A*x -b;
norm(result)

A1 = A;
v = 1:size(A1,2);
c = nchoosek(v,2);
for i = 1:size(A1,2)
    NewCol = abs(A1(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A1 = [A1 NewCol];
end
[Q, R] = custom_HQR(A1);
fact_A = Q* R;
x = R\(Q\b);
result = fact_A*x -b;
norm(result)

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
[Q, R] = custom_HQR(A2);
fact_A = Q* R;
x = R\(Q\b);
result = fact_A*x -b;
norm(result)

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
[Q, R] = custom_HQR(A3);
fact_A = Q* R;
x = R\(Q\b);
result = fact_A*x -b;
norm(result)

% title('Conjugate Gradient')
% xlabel('Iteration')
% ylabel('Loss')
% legend('Original A', 'A & Squares Cols', 'A & Colwise Mul', 'A & Squares & Colwise')