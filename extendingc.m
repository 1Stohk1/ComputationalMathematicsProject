function [A] = extendingc(A, type)
    if nargin < 2
        type = 0;
    end
    range_c= 1:size(A, 2);
    combo = nchoosek(range_c, 2);
if type == 0 || type == 1
    for i = 1:size(A, 2)
        squa_col = A(:, i).^2;
        A = [A squa_col];
%         log_col = log(A(:, i));
%         A = [A log_col];
    end
end
if type == 0 || type == 2
    for i = combo:size(combo, 1)
        j = combo(i, 1);
        k = combo(i, 2);
        mul_col = A(:, j).*A(:, k);
        A = [A mul_col];
    end
end