function[Qb, R1, x] = custom_opt_HQR(A, b)

[m, n] = size(A);
Qb = b;
for k = 1:n
    z = A(k:m, k);
    v = [-sign(z(1))*norm(z) - z(1); -z(2:end) ];
    v = v/ sqrt(v'*v);
    for j =1:n
        A(k:m,j) = A(k:m,j) - v*(2*(v'*A(k:m,j)) );
    end
    for j = 1:size(b,2)
        Qb(k:m,j) = Qb(k:m,j) - v*(2*(v'*Qb(k:m,j)) );
    end
end
R = triu(A);
R1 = R(1:n,:);
x = R1\(eye(size(R1,2), size(Qb,1))*Qb);
