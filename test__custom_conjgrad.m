temp = csvread('ML-CUP22-TR.csv', 8);
A = temp(:,2:10);
b = temp(:,11:12);
v = 1:size(A,2);
c = nchoosek(v,2);
for i = 1:size(A,2)
    NewCol = abs(A(:,i)).^2;
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k));
    NewCol = NewCol/norm(NewCol);
    A = [A NewCol];
end
clear NewCol i temp v c

b = A'*b;
A = A'*A;

init_x = zeros(size(A, 1), size(b, 2));
initial = A*init_x-b;
norm(initial)

[c, v] = custom_conjgrad(A, b(:, 1));
[d, v] = custom_conjgrad(A, b(:, 2));

result = A*[c d]-b;
norm(result)