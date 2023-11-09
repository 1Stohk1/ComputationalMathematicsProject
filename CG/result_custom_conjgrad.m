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
clear NewCol i temp v c
symm_b = A'*b;
symm_A = A'*A;

tic
[x, residual] = custom_conjgrad(symm_A, symm_b);
toc
norm(symm_A*x -symm_b)

tic
[ma,na]=size(symm_A);
[mb,nb]=size(symm_b);
afun=@(x)  reshape(symm_A*reshape(x,na,[]),[],1);
X=pcg(afun,symm_b(:));
X=reshape(X,na,nb);
toc
norm(symm_A*X -symm_b)

tic
x = symm_A\symm_b;
toc
norm(symm_A*x -symm_b)
