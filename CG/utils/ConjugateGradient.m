temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);

[new_A, new_b] = data_prep(3);

fc1 = @(x) our_case( new_A , new_b(:,1) , x );
[res1, status] = NCG(fc1, zeros(54, 1), 3, 0, 1e-6, 1, 1000);


fc2 = @(x) our_case( new_A , new_b(:,2) , x );
[res2, status] = NCG(fc2, zeros(54, 1), 3, 0, 1e-6, 1, 1000);
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






