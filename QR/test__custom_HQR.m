addpath(genpath([fileparts(pwd), filesep]));

temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);
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

[Q, R] = custom_HQR(A);
fact_A = Q* R;

if norm(A-fact_A, Inf) / norm(A, Inf) < 1e-10
    disp('The factorized matrix is equal to the first one')
else
    disp('The factorized matrix is NOT equal to the first one')
    [~,ia] = setdiff(A,Q * R);
end

x = R\(Q\b);
result = A*x -b;
norm(result)

% plot the eig of original A and R
% semilogy(abs(diag(R)),"-o")
% hold on
% semilogy(svd(A),"r-o")
% legend("Diagonal of R","Singular Values of A")

[Q,R] = qr(A);
x = R\(Q\b);
norm(A*x-b)

% The result is equal to the library!
