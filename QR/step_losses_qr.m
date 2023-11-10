addpath(genpath([fileparts(pwd), filesep]));

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);
size(A)
v = 1:size(A,2);
c = nchoosek(v,2);
for i = 1:size(A,2)
    NewCol = abs(A(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
clear NewCol i temp v c j k
tic
[m, n] = size(A);
atIter = 0;
residuals = zeros(1, n); 
for k = 1:n
    atIter = atIter + 1;
    [Q, R] = custom_HQR(A(:, 1:k));
    x = R\(Q\b);
    result = A(:, 1:k)*x -b;
    residuals(atIter) = norm(result);
end
toc
plot(residuals, 'LineWidth', 2);

title('Householder QR')
xlabel('Dimension of A')
ylabel('Loss')
% legend('Original A', 'A & Squares Cols', 'A & Colwise Mul', 'A & Squares & Colwise')
