% test__NCG
addpath(genpath([fileparts(pwd), filesep]));
temp = csvread('ML-CUP22-TR.csv', 8);

A = temp(:,2:10);
b = temp(:,11:12);

disp('ORIGINAL MATRIX')
symm_b = A'*b;
symm_A = A'*A;

% X=NCG(@(x)  genericquad(symm_A, symm_b, x), []);


% afun=@(x) symm_A*x-symm_b;
% xStar = symm_A \ symm_b;
% afun(xStar)

aFun = @(x)  genericquad(symm_A, symm_b, x);
aFun(xStar)


TF{ 6 } = @rosenbrock;

funny = @(x) rosenbrock(x);
funny([1 3])
% plotQ(symm_A, symm_b, [-100000 100000])

X=NCG(@(x) rosenbrock(x), [], 2);

function plotQ( Q , q , range)

f = @(x,y) 0.5 * [ x y ] * Q * [ x ; y ] + q' * [ x ; y ];

warning( 'off' , 'all' );

fcontour( f , range , 'LineColor' , 'k' , 'LineWidth' , 1 );

warning( 'on' , 'all' );

end

function [ v , varargout ] = genericquad( Q , q , x )

if isempty( x )  % informative call
   if min( eig( Q ) ) > 1e-14
      xStar = Q \ q;
      v = Q * xStar - q;
   else
      v = - Inf;
   end
   if nargout > 1
      varargout{ 1 } = [ 0 ; 0 ];
   end
else
   v = Q * x - q;  % f(x)
   if nargout > 1
      varargout{ 1 } = x;  % \nabla f(x)
      if nargout > 2
         varargout{ 2 } = x;       % \nabla^2 f(x)
      end
   end
end
end  % genericquad


function [ v , varargout ] = rosenbrock( x )
% rosenbrock's valley-shaped function
% syms x y
% f = @(x, y) 100 * ( y - x^2 )^2 + ( x - 1 )^2
%
% diff( f , x )
% 2 * x - 400 * x * ( - x^2 + y ) - 2
%
% diff( f , y )
% - 200 * x^2 + 200 * y
%
% diff( f , x , 2 )
% 1200 * x^2 - 400 * y + 2
%
% diff( f , y , 2 )
% 200
%
% diff( f , x , y )
% -400 * x

if isempty( x )  % informative call
   v = 0;
   if nargout > 1
      varargout{ 1 } = [ -1 ; 1 ];
   end
else
   v = 100 * ( x( 2 ) - x( 1 )^2 )^2 + ( x( 1 ) - 1 )^2;  % f(x)
   if nargout > 1
      g = zeros( 2 , 1 );
      g( 1 ) = 2 * x( 1 ) - 400* x( 1 ) * ( x( 2 ) - x( 1 )^2 ) - 2;
      g( 2 ) = - 200 * x( 1 )^2 + 200 * x( 2 );
      
      varargout{ 1 } = g;  % \nabla f(x)
      if nargout > 2
         H = zeros( 2 , 2 );
         H( 1 , 1 ) = 1200 * x( 1 )^2 - 400 * x( 2 ) + 2;
         H( 2 , 2 ) = 200;
         H( 2 , 1 ) = -400 * x( 1 );
         H( 1 , 2 ) = H( 2 , 1 );
         varargout{ 2 } = H;       % \nabla^2 f(x)
      end
   end
end
end  % rosenbrock