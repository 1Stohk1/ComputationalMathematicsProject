function[x,  status] = alternative_conjgrad(A, b, x, tol)

    functionRows = size(A, 1);
    functionCols = size(b, 2);
    
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
    temp_z = zeros(functionRows, functionCols); 
    temp_s = zeros(1, functionCols); 

    for i = 1:functionCols
        residual(:, i) = b(:, i) - A * x(:, i);
        direction(:, i) = -residual(:, i);
        temp_z(:, i)  = A*direction(:, i);
        temp_s(i) = direction(:, i)'*temp_z(:, i);
        alpha(i) = (residual(:, i)'*direction(:, i))/temp_s(i);
        x(:, i) = x(:, i) + alpha(i)* direction(:, i);
    end

%    Starting the loop
    while norm(residual) > tol  && atIter < maxIters 

        for i = 1:functionCols
            residual(:, i) = residual(:, i) - alpha(i)*temp_z(:, i);
            beta(i) = (residual(:, i)'*temp_z(:, i))/temp_s(i);
            direction(:, i) = -residual(:, i) + beta(i)*direction(:, i);
            
            temp_z(:, i)  = A*direction(:, i);
            temp_s(i) = direction(:, i)'*temp_z(:, i);

            alpha(i) = (residual(:, i)'*direction(:, i))/temp_s(i);
            x(:, i) = x(:, i) + alpha(i)* direction(:, i);
        end
        
        atIter   = atIter + 1;
        resVals(atIter) = norm(A*x-b);
        
    end
    %      plot(resVals)
    loglog(resVals)

    status = 0;

