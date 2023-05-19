temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

v = 1:size(A,2);
c = nchoosek(v,2);

disp('Condition number of initial A')
cond(A)

for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    NewCol1 = A(:,j).^2;
    NewCol1 = NewCol1/norm(NewCol1);
    A = [A NewCol NewCol1];
end

disp('Condition number of A after Column additioning')
cond(A)

sim_A = A'*A;
D = diag(sim_A);

C = sqrt(inv(diag(D)))*sim_A*sqrt(inv(diag(D)));
alpha = min(eig(C));
%alpha = 1e-16

new_A = sim_A+alpha*D;
%new_A = sim_A;

new_b = (A'*b)+alpha*D;
%new_b = A'*b;

disp('Condition number of final A')
cond(new_A)

try all(eig(new_A) > 0)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

x1 = pcg(new_A,new_b(:,1));
x2 = pcg(new_A,new_b(:,2));
x = [x1 x2];

residual = new_A*x - new_b;

sqrt(norm(residual))
if sqrt(norm(residual))<1000
    disp('Bella roba')
end

