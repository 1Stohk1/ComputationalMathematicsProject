function[x,  resVals, atIter, exTime] = custom_conjgrad(A, b, x, tol, maxIters, metric, A_original, b_original)

% Application of the Conjugate Gradient to the inserted system formed by A and b. Can be 
% applied to systems with more than 1 column goal vector b (the algorithm will iterate for 
% each column in b).
% 
% -- Input Arguments:
% - x: starting point of the solution vector;
% - tol: tolerance that will break the main loop;
% - maxIters: maximum number of iteration performed by the algorithm;
% - metric: metric to be followed in order to break the loop from the main loop.
% 
% -- Output Arguments:
% - x: solution vector;
% - resVals: vector with each residual computed;
% - atIter: scalar that represent the number of iterations performed;
% - exTime: execution time for each iteration, the item at end will represent the sum of all
% the precedent iterations.

    [functionRows, functionCols] = size(b);
    
    %     Validation of the input args
    if nargin<3
        x = zeros(functionRows, functionCols);
    end
    if nargin<4
        tol = 1e-8;
    end
    if nargin<5
        maxIters = functionRows;
    end    
    if nargin<6
        metric = 0;
    end
    if nargin<7
        A_original = A;
    end
    if nargin<8
        b_original = b;
    end

    if  ~ isequal(size(x), size(b))
        warning(['The input x must have the same dimension of b, size(x): (%d, %d).' ...
            '\nTherefore default zeros(size(b)) will be used.'], size(x))
        x = zeros(functionRows, functionCols);
    end
    if maxIters <1
        warning(['At least 1 iteration must be performed, inserted: %d.' ...
            '\nTherefor default iteration equal to numCols will be used.'], maxIters)
        maxIters = functionCols;
    end
    if metric <0 && modality >1
        warning(['There are only 2 metrics that can be watched resVals in modality 0 (default) and diff_star in modality 1, inserted: %d.' ...
            '\nTherefore default will be used'], metric)
        metric = 0;
    end
    if metric == 1
        x_star = A_original\b_original;
    end

    
    %     Initialize the variables
    tic;
    atIter = 1;
    pre = atIter;
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
    
    exTime(atIter) = toc;
    
%    Starting the loop
    while resVals(pre) > tol  && atIter <= maxIters
        tic;

        for i = 1:functionCols
            new_residual(:, i) = residual(:, i) - alpha(i)*A*direction(:, i);
            beta(i) = (new_residual(:, i)'*new_residual(:, i))/(residual(:, i)'*residual(:, i));
            
            residual(:, i) = new_residual(:, i);
            
            direction(:, i) = new_residual(:, i) + beta(i)*direction(:, i);
            alpha(i) = (new_residual(:, i)'*new_residual(:, i))/(direction(:, i)'*A*direction(:, i));

            x(:, i) = x(:, i) + alpha(i)* direction(:, i);
        end
        
        if metric == 0
            resVals(atIter) = norm(A_original*x-b_original)/norm(b_original);
        else 
            resVals(atIter) = norm(x-x_star)/norm(x_star);
        end

        exTime(atIter) = exTime(atIter)+toc;
        pre = atIter;
        atIter = atIter + 1;
        exTime(atIter) = exTime(pre);
    end

