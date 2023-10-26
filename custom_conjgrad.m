function[x,  status] = custom_conjgrad(A, b, x, tol)

    [functionRows, functionCols] = size(b);
    
    %     Control the input args
    if nargin<3
        x = zeros(functionRows, functionCols);
    end
    if nargin<4
        tol = 1e-6;
    end

    %     Initialize the variables
    atIter = 0;
    maxIters = functionRows;
    resVals = zeros(1, maxIters); 
   
    direction = zeros(functionRows, functionCols); 
    residual = zeros(functionRows, functionCols); 
    new_residual = zeros(functionRows, functionCols); 
    alpha = zeros(1, functionCols); 
    beta = zeros(1, functionCols); 
    
    for i = 1:functionCols
        direction(:, i) = b(:, i) - A * x(:, i);
        residual(:, i) = direction(:, i);
        alpha(i) = (residual(:, i)'*residual(:, i))/(direction(:, i)'*A*direction(:, i));
        x(:, i) = x(:, i) + alpha(i)* direction(:, i);
    end

%    Starting the loop
    while norm(residual) > tol  && atIter < maxIters 

        for i = 1:functionCols
            new_residual(:, i) = residual(:, i) - alpha(i)*A*direction(:, i);
            beta(i) = (new_residual(:, i)'*new_residual(:, i))/(residual(:, i)'*residual(:, i));
            
            residual(:, i) = new_residual(:, i);
            
            direction(:, i) = new_residual(:, i) + beta(i)*direction(:, i);
            alpha(i) = (new_residual(:, i)'*new_residual(:, i))/(direction(:, i)'*A*direction(:, i));

            x(:, i) = x(:, i) + alpha(i)* direction(:, i);
        end
        
        atIter   = atIter + 1;
        resVals(atIter) = norm(A*x-b);
        
    end
    %      plot(resVals)
    loglog(resVals)

    status = 0;

