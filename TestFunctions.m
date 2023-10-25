function TF = TestFunctions()

%function TF = TestFunctions()
%
% Produces a cell array of function handlers, useful to test unconstrained
% optimization algorithms.
%
% Each function in the array has the following interface:
%
%   [ v , varargout ] = f( x )
%
% Input:
%
% - x is either a [ n x 1 ] real (column) vector denoting the input of
%   f(), or [] (empty).
%
% Output:
%
% - v (real, scalar): if x == [] this is the best known lower bound on
%   the unconstrained global optimum of f(); it can be -Inf if either f()
%   is not bounded below, or no such information is available. If x ~= []
%   then v = f(x).
%
% - g (real, [ n x 1 ] real vector) is the first optional argument. This
%   also depends on x. if x == [] this is the standard starting point of an
%   optimization algorithm, otherwise it is the gradient of f() at x, or a
%   subgradient if f() is not differentiable at x.
%
% - H (real, [ n x n ] real matrix) is the first optional argument. This
%   must only be specified if x ~= [], and it is the Hessian of f() at x.
%   If no such information is available, the function throws error.
%
% The current list of functions is the following:
%
%  1 Standard 2x2 PSD quadratic function with nicely conditioned Hessian.
%
%  2 Standard 2x2 PSD quadratic function with less nicely conditioned
%    Hessian.
%
%  3 Standard 2x2 PSD quadratic function with Hessian having one zero
%    eigenvalue.
%
%  4 Standard 2x2 quadratic function with indefinite Hessian (one positive
%    and one negative eigenvalue)
%
%  5 Standard 2x2 quadratic function with "very elongated" Hessian (a 
%    very small positive minimum eigenvalue, the other much larger)
%
%  6 the 2-dim Rosenbrock function
%
%  7 the "six-hump camel" function
%
%  8 the Ackley function
%
%  9 a 2-dim nondifferentiable function coming from Lasso regularization
%
% 10 a 76-dim (nonconvex, differentiable) function coming from a fitting
%    problem with ( X , y ) both [ 288 , 1 ] (i.e., a fitting with only
%    one feature) using a "rough" NN with 1 input, 1 output, 3 hidden
%    layers of 5 nodes each, and tanh activation function
%
% 11 same as 10 plus a 1e-4 || x ||^2 / 2 ridge stabilising term
%
%{
 =======================================
 Author: Antonio Frangioni
 Date: 08-11-18
 Version 1.01
 Copyright Antonio Frangioni
 =======================================
%}

TF = cell( 10 , 1 );
TF{ 1 } = @(x) genericquad( [ 6 -2 ; -2 6 ] , [ 10 ; 5 ] , x );
% eigenvalues: 4, 8
TF{ 2 } = @(x) genericquad( [ 5 -3 ; -3 5 ] , [ 10 ; 5 ] , x );
% eigenvalues: 2, 8
TF{ 3 } = @(x) genericquad( [ 4 -4 ; -4 4 ] , [ 10 ; 5 ] , x );
% eigenvalues: 0, 8
TF{ 4 } = @(x) genericquad( [ 3 -5 ; -5 3 ] , [ 10 ; 5 ] , x );
% eigenvalues: -2, 8
TF{ 5 } = @(x) genericquad( [ 101 -99 ; -99 101 ] , [ 10 ; 5 ] , x );
% eigenvalues: 2, 200
% HBG: alpha = 0.0165 , beta = 0.678
TF{ 6 } = @rosenbrock;
TF{ 7 } = @sixhumpcamel;
TF{ 8 } = @ackley;
TF{ 9 } = @lasso;
TF{ 10 } = @myNN;
TF{ 11 } = @myNN2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = genericquad( Q , q , x )
 % generic quadratic function f(x) = x' * Q * x / 2 + q' * x

if isempty( x )  % informative call
   if min( eig( Q ) ) > 1e-14
      xStar = Q \ -q;
      v = 0.5 * xStar' * Q * xStar + q' * xStar;
   else
      v = - Inf;
   end
   if nargout > 1
      varargout{ 1 } = [ 0 ; 0 ];
   end
else
   if ~ isequal( size( x ) , [ 2 1 ] )
      error( 'genericquad: x is of wrong size' );
   end
   v = 0.5 * x' * Q * x + q' * x;  % f(x)
   if nargout > 1
      varargout{ 1 } = Q * x + q;  % \nabla f(x)
      if nargout > 2
         varargout{ 2 } = Q;       % \nabla^2 f(x)
      end
   end
end
end  % genericquad

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = sixhumpcamel( x )
% six-hump-camel valley-shaped function
% syms x y
% f = @(x, y) ( 4 - 2.1 * x^2 + x^4 / 3 ) * x^2 + x * y + 4 * ( y^2 - 1 ) *
% y^2
%
% diff( f , x )
% 2 * x^5 - ( 42 * x^3 ) / 5 + 8 * x + y
%
% diff( f , y )
% 16 * y^3 - 8 * y + x
%
% diff( f , x , 2 )
% 10 * x^4 - ( 126 * x^2 ) / 5 + 8
%
% diff( f , y , 2 )
% 48 * y^2 - 8
%
% diff( f , x , y )
% 1

if isempty( x )  % informative call
   v = -1.03162845349;
   if nargout > 1
      varargout{ 1 } = [ 0 ; 0 ];
   end
else
   v = ( 4 - 2.1 * x( 1 )^2 + x( 1 )^4 / 3 ) * x( 1 )^2 + ...
       x( 1 ) * x( 2 ) + 4 * ( x( 2 )^2 - 1 ) * x( 2 )^2;  % f(x)
   if nargout > 1
      g = zeros( 2 , 1 );
      g( 1 ) = 2 * x( 1 )^5 - ( 42 * x( 1 )^3 ) / 5 + 8 * x( 1 ) + x( 2 );
      g( 2 ) = 16 * x( 2 )^3 - 8 * x( 2 ) + x( 1 );
      
      varargout{ 1 } = g;  % \nabla f(x)
      if nargout > 2
         H = zeros( 2 , 2 );
         H( 1 , 1 ) = 10 * x( 1 )^4 - ( 126 * x( 1 )^2 ) / 5 + 8;
         H( 2 , 2 ) = 48 * x( 2 )^2 - 8;
         H( 2 , 1 ) = 1;
         H( 1 , 2 ) = H( 2 , 1 );
         varargout{ 2 } = H;       % \nabla^2 f(x)
      end
   end
end
end  % sixhumpcamel

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = ackley( xx )

% syms x y
% f = @(x, y) - 20 * exp( - 0.2 * sqrt( ( x^2 + y^2 ) / 2 ) ) ...
%             - exp( ( cos( 2 * pi * x ) + cos( 2 * pi * y ) ) / 2 ) ...
%             + 20 + exp(1)
%

ManuallyComputedfGH = 0;

if isempty( xx )  % informative call
   v = 0;
   if nargout > 1
      varargout{ 1 } = [ 2 ; 2 ];
   end
else
   if ~ isequal( size( xx ) , [ 2 1 ] )
      error( 'ackley: x is of wrong size' );
   end

   if ManuallyComputedfGH

% diff( f , x )
% pi*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*sin(2*pi*x) +
% (2*x*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(1/2)
%
% diff( f , y )
% pi*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*sin(2*pi*y) +
% (2*y*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(1/2)
%
% diff( f , x , 2 )
% 
% (2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(1/2) +
% 2*pi^2*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*cos(2*pi*x) -
% (x^2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(5*(x^2/2 + y^2/2)) -
% (x^2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(3/2) -
% pi^2*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*sin(2*pi*x)^2
%
% diff( f , y , 2 )
% (2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(1/2) +
% 2*pi^2*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*cos(2*pi*y) -
% (y^2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(5*(x^2/2 + y^2/2)) -
% (y^2*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(3/2) -
% pi^2*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*sin(2*pi*y)^2
%
% diff( f , x , y)
% - (x*y*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(5*(x^2/2 + y^2/2)) -
% (x*y*exp(-(x^2/2 + y^2/2)^(1/2)/5))/(x^2/2 + y^2/2)^(3/2) -
% pi^2*exp(cos(2*pi*x)/2 + cos(2*pi*y)/2)*sin(2*pi*x)*sin(2*pi*y)

      x = xx( 1 );
      y = xx( 2 );
      sqn2 = ( x^2 + y^2 ) / 2;
      cosx = cos( 2 * pi * x );
      cosy = cos( 2 * pi * y );
      comp1 = exp( - (sqn2)^(1/2) / 5 );
      comp2 = exp( ( cosx + cosy ) / 2 );

      v = - 20 * comp1 - comp2 + 20 + exp( 1 );  

      if nargout > 1
         sinx = sin( 2 * pi * x );
         siny = sin( 2 * pi * y );

         g = zeros( 2 , 1 );
         g( 1 ) = pi * comp2 * sinx + 2 * x * comp1 / (sqn2)^(1/2);
         
         g( 2 ) = pi * comp2 * siny + 2 * y * comp1 / (sqn2)^(1/2);
         
         varargout{ 1 } = g;  % \nabla f(x)
         if nargout > 2
            H = zeros( 2 , 2 );

            H( 1 , 1 ) = (2*comp1)/(sqn2)^(1/2) + 2*pi^2*comp2*cosx     ...
                       - (x^2*comp1)/(5*sqn2) - (x^2*comp1)/(sqn2)^(3/2)...
                       - pi^2*comp2*sinx^2;

            H( 2 , 2 ) = (2*comp1)/(sqn2)^(1/2) + 2*pi^2*comp2*cosy     ...
                       - (y^2*comp1)/(5*sqn2) - (y^2*comp1)/(sqn2)^(3/2)...
                       - pi^2*comp2*siny^2;

            H( 1 , 2 ) = - (x*y*comp1)/(5*(sqn2))                       ...
                         - (x*y*comp1)/(sqn2)^(3/2)                     ...
                         - pi^2*comp2*sinx*siny;

            H( 2 , 1 ) = H( 1 , 2 );
            varargout{ 2 } = H;       % \nabla^2 f(x)
         end
      end
   else
      if nargout > 2
         [ H , g , v ] = ackley_Hes( xx );
         varargout{ 2 } = H;
         varargout{ 1 } = g';
      elseif nargout > 1
         [ g , v ] = ackley_Grd( xx );
         varargout{ 1 } = g';
      else
         v = - 20 * exp( - ( ( xx( 1 )^2 + xx( 2 )^2 ) / 2 )^(1/2) / 5 )...
             - exp( cos( 2 * pi * xx( 1 ) ) / 2 + ...
                    cos( 2 * pi * xx( 2 ) ) / 2 ) + 20 + exp( 1 );  
      end
   end
end
end  % ackley

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = lasso( x )
% nondifferentiable lasso example:
%
% f( x , y ) = || 3 * x + 2 * y - 2 ||_2^2 + 10 ( | x | + | y | )

if isempty( x )  % informative call
   v = ( 2 - 1/3 )^2 + 10/9;  % optimal solution [ 1/9 , 0 ]
   if nargout > 1
      varargout{ 1 } = [ 0 ; 0 ];
   end
else
   v = ( 3 * x( 1 ) + 2 * x( 2 ) - 2 )^2 + ...
                           10 * ( abs( x( 1 ) ) + abs( x( 2 ) ) );  % f(x)  
   if nargout > 1
      g = zeros( 2 , 1 );
      g( 1 ) = 18 * x( 1 ) + 12 * x( 2 ) - 12 + 10 * sign( x( 1 ) );
      g( 2 ) = 12 * x( 1 ) + 8 * x( 2 ) - 8 + 10 * sign( x( 2 ) );
      
      varargout{ 1 } = g;  % \nabla f(x)
      if nargout > 2
         error( 'lasso: Hessian not available' );
      end
   end
end
end  % lasso

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = myNN( x )
% 1 x 5 x 5 x 5 x 1 = 76 w NN for solving a 1D fitting problem

if isempty( x )  % informative call
   v = - Inf;  % optimal value unknown (although 0 may perhaps be good)
   if nargout > 1
      % Xavier initialization: uniform random in [ - A , A ] with
      % A = \sqrt{6} / \sqrt{n + m}, with n and m the input and output
      % layers. in our case n + m is either 6 or 10, so we take A = 1
      %
      % note that starting point is random, so each run will be different
      % (unless an explicit starting point is provided); if stability is
      % neeed, the seed of the generator has to be set externally
      varargout{ 1 } = 2 * rand( 76 , 1 ) - 1;
   end
else
   v = testNN( x );  % f(x)  
   if nargout > 1     
      varargout{ 1 } = testNN_Jac( x )';     % \nabla f( x )
      if nargout > 2
         varargout{ 2 } = testNN_Hes( x )';  % \nabla^2 f( x )
      end
   end
end
end  % myNN

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function [ v , varargout ] = myNN2( x )
% 1 x 5 x 5 x 5 x 1 = 76 w NN for solving a 1D fitting problem
% plus ridge stabilization \lambda || x ||^2 / 2

lambda = 1e+2;

if isempty( x )  % informative call
   v = - Inf;  % optimal value unknown (although 0 may perhaps be good)
   if nargout > 1
      % Xavier initialization: uniform random in [ - A , A ] with
      % A = \sqrt{6} / \sqrt{n + m}, with n and m the input and output
      % layers. in our case n + m is either 6 or 10, so we take A = 1
      %
      % note that starting point is random, so each run will be different
      % (unless an explicit starting point is provided); if stability is
      % neeed, the seed of the generator has to be set externally
      varargout{ 1 } = 2 * rand( 76 , 1 ) - 1;
   end
else
   v = testNN( x ) + lambda * x' * x / 2;  % f(x)  
   if nargout > 1     
      varargout{ 1 } = testNN_Jac( x )' + lambda * x;  % \nabla f( x )
      if nargout > 2
         varargout{ 2 } = testNN_Hes( x )' + lambda * eye( 76 );
                                                       % \nabla^2 f( x )
      end
   end
end
end  % myNN2

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end