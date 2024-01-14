function [A, b] = data_prep(method, type)
    if nargin <1
        warning(['You should choose at least one from:' ...
        '\n1) Tall thin (0);'...
        '\n2) Symmetric (1).'])
        method = 0;
    end
    if nargin <2
        type = 0;
    elseif type < -1 || type > 2
        warning(['The extension combinations are: ', ...
        '\n1) Not applying extensions (-1); \n2) Performing all the extensions (0);', ...
        '\n3) Appending only the square of the columns (1);', ...
        '\n4) Appending only the multiplication between columns (2).', ...
        '\nInserted: %.f, hence default 0 will be used'], type);
        type = 0;
    end

    temp = csvread('ML-CUP22-TR.csv', 8);
    A = temp(:, 2:10);
    b = temp(:, 11);
    if ne(type, -1)
        A = extendingc(A, type);
    end
    A = normc(A);
    if method == 1
        b = A'*b;
        A = A'*A;
    end
end
