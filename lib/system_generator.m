function[A, b] = system_generator(matrixColumn, matrixRow, b_col)
    rng(42)
    while true
      A = rand(matrixRow, matrixColumn);
      if rank(A) == matrixColumn; break; end    %will be true nearly all the time
    end
    b = rand(matrixRow, b_col);
end