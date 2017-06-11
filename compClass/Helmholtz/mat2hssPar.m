function [ Q1, Q2, L, DL, U, B, W, V, S, nflops, storage ] = mat2hssPar( A, hssTT, hssTL, rg, tol, SchurInd )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  HSS parallel compression, factorization, Schur Complement %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nflops = 0;

% 1. HSS compression
[ D, U, R, B, W, V, nflops1 ] = mat2hssParcmp( A, hssTT, hssTL, rg, tol );
nflops = nflops + nflops1; clear nflops1;


% 2. HSS factorization
[ Q1, Q2, L, DL, U, B, W, V, S, nflops1, storage ] = mat2hssParfac( D, U, R, B, W, V, hssTT, hssTL, SchurInd );
nflops = nflops + nflops1; clear nflops1;


return;
end