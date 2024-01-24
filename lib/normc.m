function [A] = normc(A)
    for i = 1:size(A,2)
        col = A(:,i)/norm(A(:,i));
        A(:,i) =col;
    end
end