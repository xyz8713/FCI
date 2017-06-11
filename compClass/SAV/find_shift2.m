function [t] = find_shift2(A)
%% UNUSED for now
%% finds two shifts so that
%% A - t_1 I has only Real(lamda)>0 eigenvalues 
%% A - t_2 I has only real(lambda)<0 eigenvalues 
%% A may be complex
%% sum of abs values in each row
n = size(A,1);
d = abs(A)*ones(n,1) - abs(diag(A));
n = size(A,1);
%% first part
d1 = min(real(diag(A))-d);
d2 = max(real(diag(A))-d);
t = [d1, d2];
disp(t);
