function [t] = find_shift(A)
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
t = 0;

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
    alp = -bet + s*del;
    if (abs(alp) > abs(t)) 
        t = alp;
    end
end
%% t = t*0.25i;
   t = t*0.75i;
   disp(t);
    
