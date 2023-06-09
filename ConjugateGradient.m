temp = csvread("Documenti/Matlab Scripts/ML-CUP22-TR.csv", 8);
A = temp(:,2:10);
b = temp(:,11);

v = 1:size(A,2);
c = nchoosek(v,2);

disp('Condition number of initial A')
cond(A)

for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end

disp('Condition number of A after Column additioning')
cond(A)

new_A = A'*A;
new_b = A'*b;

disp('Condition number of final A')
cond(new_A)

try all(eig(new_A) > 0)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

x1 = pcg(new_A,new_b(:,1));
%x2 = pcg(new_A,new_b(:,2));
%x = [x1 x2];

residual = new_A*x1 - new_b;

sqrt(norm(residual))
if sqrt(norm(residual))<1000
    disp('Bella roba')
end
