function[Q, R1, x] = custom_opt_HQR(A, b)

[m, n] = size(A);
Q = b;
for k = 1:n
    z = A(k:m, k);
    v = [-sign(z(1))*norm(z) - z(1); -z(2:end) ];
    v = v/ sqrt(v'*v);
    for j =1:n
        A(k:m,j) = A(k:m,j) - v*(2*(v'*A(k:m,j)) );
    end
    for j = 1:size(b,2)
        Q(k:m,j) = Q(k:m,j) - v*(2*(v'*Q(k:m,j)) );
    end
end
R = triu(A);
R1 = R(1:n,:);
x = R1\(eye(size(R1,2), size(Q,1))*Q);
