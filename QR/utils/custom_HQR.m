function[Q, R] = custom_HQR(A)

[m, n] = size(A);
Q = eye(m);
for k = 1:n
%     Take the column vector that we want to reflect 
%     we care only from the k till the end rows cause all the precedent are
%     already cancelled, and we took the k column
    z = A(k:m, k);
%     we want to preserve the first element so we take the opposite of the
%     first element, multiplied by its length, such a way we will not delete him, 
%     and attach the rest of the vector.
    v = [-sign(z(1))*norm(z) - z(1); -z(2:end) ];
%     This will be our vector in which we build the mirror, such a way we can
%     cancel everything (cause we will subtract the parallel stuff to this, i.e. 
%     everything cause the mirror is just doing this
    v = v/ norm(v);
%     The length of w (vector to the mirror) is the length of x * cosine
%     theta, and this length can be computed with the length of negative x
%     and v (the vector in which we built the orthogonal mirror), the
%     direction will be the same of v so the vector w is equal to its
%     length and direction multiplied so 
%     -x'*v/sqrt(v'*v)*v/sqrt(v'*v)  = -v*(x'*v)/(v'*v)
    for j =1:n
        A(k:m,j) = A(k:m,j) - v*(2*(v'*A(k:m,j)) );
    end
%     In order to get the reflection we will add the double length of
%     the vector pointing the mirror (one for the mirror another one for
%     the reflection, we will apply the reflection to both A and Q in order
%     to get track of the column vectors to delete the bottom triangle of A
%     and get R from A
    for j = 1:m
        Q(k:m,j) = Q(k:m,j) - v*(2*(v'*Q(k:m,j)) );
    end
end

R = triu(A);
