 function x = precLU(PRE, rhs) 
%% function rhs = precLU(PRE, rhs) 
%% arms preconditioning operation
%% PRE = struct for preconditioner
%%-------------------------------------------------
L = PRE.L;
U = PRE.U; 
P = PRE.P;
x = U \ (L \ (P*rhs)); 

