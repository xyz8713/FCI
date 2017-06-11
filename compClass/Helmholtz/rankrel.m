function r = rankrel(A,tol)

s = svd(full(A));
r = length(find(s >= s(1)*tol));