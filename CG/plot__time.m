addpath(genpath([fileparts(pwd), filesep]));

dim = 100;
[A, b] = system_generator(dim, dim, 1);

b = A'*b;
A = A'*A;
times = zeros(1, dim);
for j= 1:1000
    [x, res, iter, ex_time] = custom_conjgrad(A, b, zeros(size(b)), 0);
    times = times + ex_time;
end
% Plotting the matrix + colwise multiplication results
times = times*1000;
plot(times, 'r-','LineWidth',2);
xlim([0 dim])

title('Conjugate Gradient (Time)')
xlabel('Iterations')
ylabel('Ex. Time ($\mu$)', 'Interpreter', 'Latex')