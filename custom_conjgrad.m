function[x,  status] = custom_conjgrad(A, b, x, tol)

   atIter   = 0;
   maxIters = numel(b);
   resVals    = zeros(1, maxIters); 


%     Control on the number of input args
   if nargin<3
        x = zeros(size(A, 1), 1);
        size(x)
   end
   if nargin<4
       tol = 1e-6;
   end
    
    direction = b - A * x;
    size(direction)
    residual = direction;
    alpha = (residual'*residual)/(direction'*A*direction);
    size(alpha)
     x = x + alpha* direction;
    
    while norm(residual) > tol  && atIter < maxIters 
        
       new_residual = residual - alpha*A*direction;
       beta = (new_residual'*new_residual)/(residual'*residual);
       direction = new_residual + beta*direction;
       alpha = (new_residual'*new_residual)/(direction'*A*direction);
       
        x = x + alpha* direction;
        residual = new_residual;
        atIter   = atIter + 1;
        resVals(atIter) = norm(residual);
    end
%      plot(resVals)
     loglog(resVals)

    status = 0;

