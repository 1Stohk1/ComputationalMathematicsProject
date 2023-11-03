temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

v = 1:size(A,2);
c = nchoosek(v,2);

disp(['Condition number of initial A ', num2str(cond(A))])

disp(['Rank of initial A ',  num2str(rank(A)), ' sizes (' , num2str(size(A)), ')'])
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
new_A = A'*A;
new_b = A'*b;
disp(['Rank of symm A ',  num2str(rank(new_A)), ' sizes (' , num2str(size(new_A)), ')'])

if rank(new_A) ==size(new_A,1) 
    disp('Matrix is full rank.')
else
    disp('Matrix is NOT full rank')
end

disp(['Condition number of symm A ', num2str(cond(new_A))])

if all(eig(new_A) > 0)
    disp('Matrix is symmetric positive definite.')
else
    disp('Matrix is NOT symmetric positive definite')
end

x1 = conjgrad(new_A, new_b(:,1), 1e-04);
x2 = conjgrad(new_A, new_b(:,2), 1e-04);
x = [x1 x2];


residual = new_A*x - new_b;
norm(residual)

if sqrt(norm(residual))<1
    disp('Bella roba')
end

fc1 = @(x) our_case( new_A , new_b(:,1) , x );
[res1, status] = NCG(fc1, zeros(54, 1), 0, 0, 1e-6, 1, 1000);


fc2 = @(x) our_case( new_A , new_b(:,2) , x );
[res2, status] = NCG(fc, zeros(54, 1), 0, 0, 1e-6, 1, 1000);
x = [res1 res2];

residual = new_A*x - new_b;
norm(residual)


% Our function
function [ v , varargout ] = our_case( Q , q , x )
 % generic quadratic function f(x) = x' * Q * x / 2 + q' * x

if isempty( x )  % informative call
   if min( eig( Q ) ) > 1e-14
      xStar = Q \ -q;
      v = xStar' * Q * xStar + q' * xStar;
   else
      v = - Inf;
   end
   if nargout > 1
      varargout{ 1 } = [ 0 ; 0 ];
   end
else
   if ~ isequal( size( x ) , [ 54 1 ] )
      error( 'genericquad: x is of wrong size' );
   end
   v = x' * Q * x + q' * x;  % f(x)
   if nargout > 1
      varargout{ 1 } = Q * x + q;  % \nabla f(x)
      if nargout > 2
         varargout{ 2 } = Q;       % \nabla^2 f(x)
      end
   end
end
end  % genericquad






