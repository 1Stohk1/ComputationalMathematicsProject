function[Qb, R, x] = test_HQR(A, b)

% Application of the Houdeholder QR method with the solution matrix computed thanks 
% to the premultiplication of b.
% 
% -- Input Arguments
% 
% -- Output Arguments
% Qb: Orthogonal matrix premultiplied by b, this holds all the reflector vectors, that removed the components under the main diagonal of A;
% R1: The upper triangular matrix till the n-th row;
% x: Solution of the system.

[m, n] = size(A);
Qb = b;
for k = 1:n
    z = A(k:m, k);
    v = [-sign(z(1))*norm(z) - z(1); -z(2:end) ];
    v = v/ sqrt(v'*v);
    v2 = v*2;
    A(k:m,k:n) = A(k:m,k:n) - v2*(v'*A(k:m,k:n));
    Qb(k:m,:) = Qb(k:m,:) - v2*(v'*Qb(k:m,:));
end
R = triu(A(1:n,:));
x = R\(eye(size(R,2), size(Qb,1))*Qb);