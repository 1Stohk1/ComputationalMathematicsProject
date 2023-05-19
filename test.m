% Specify the file path
file = "Scaricati/Telegram Desktop/ML-CUP22-TR.csv";

% Load the entire CSV file
data = csvread(file, 8);

% Define the range of rows and columns you want to extract
startRow = 1;
endRow = size(data, 1);
startCol = 2;
endCol = size(data, 2) - 2;

% Extract the desired rows and columns
A = data(startRow:endRow, startCol:endCol);
B = data(startRow:endRow, endCol+1:size(data, 2));

v = 1:8;
c = nchoosek(v,2);

for i = c:size(c,1)
    j = c(i, 1);
    k = c(i, 2);
    NewCol = abs(A(:,j).*A(:,k)*100);
    A = [A NewCol];
end

newA = A'*A;

all(eig(newA) > 0)


all(eig(newA) < 1000)

max(eig(newA))
