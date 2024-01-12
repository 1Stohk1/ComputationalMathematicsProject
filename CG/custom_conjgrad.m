function[x,  resVals, atIter, exTime] = custom_conjgrad(A, b, x, tol, maxIters)

    [functionRows, functionCols] = size(b);
    
    %     Control the input args
    if nargin<3
        x = zeros(functionRows, functionCols);
    end
    if nargin<4
        tol = 1e-8;
    end
    if nargin<5
        maxIters = functionRows;
    end    
    
    %     Initialize the variables
    tic;
    atIter = 1;
    pre=atIter;
    resVals = zeros(1, maxIters) + 1; 
    exTime = zeros(1, maxIters); 
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
    
    exTime(1) = toc;
    
%    Starting the loop
    while resVals(pre) > tol  && atIter < maxIters 
        tic;

        for i = 1:functionCols
            new_residual(:, i) = residual(:, i) - alpha(i)*A*direction(:, i);
            beta(i) = (new_residual(:, i)'*new_residual(:, i))/(residual(:, i)'*residual(:, i));
            
            residual(:, i) = new_residual(:, i);
            
            direction(:, i) = new_residual(:, i) + beta(i)*direction(:, i);
            alpha(i) = (new_residual(:, i)'*new_residual(:, i))/(direction(:, i)'*A*direction(:, i));

            x(:, i) = x(:, i) + alpha(i)* direction(:, i);
        end
        
        resVals(atIter) = norm(A*x-b)/norm(b);
        exTime(atIter) = exTime(atIter)+toc;
        pre=atIter;
        atIter   = atIter + 1;
        exTime(atIter) = exTime(atIter-1);
    end

