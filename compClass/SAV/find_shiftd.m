function [d] = find_shiftd(A)
%% finds a shift to make A diagonally dominant
%% A is complex
%% sum of abs values in each row
n = size(A,1);
d = abs(A)*ones(n,1);
n = size(A,1);
%% See article: 
%% Daniel Ossei Kuffuor, Y. Saad. 
%% Preconditioning Helmholtz linear systems
%% tech. report. 2009 

d = d - abs(diag(A));

for i=1:n
    z = A(i,i);
    g = d(i);
    eta = real(z);
    if (abs(eta)>g) 
        continue;
    end   
    bet = imag(z) ;
     
    del = sqrt(g^2 - eta^2);
    s = 1;
    if (bet<0), s = -1.0;,  end;
    d(i) = -bet + s*del;
end   
   d = d .* 0.25i;

    
