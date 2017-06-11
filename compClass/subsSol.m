     function y = subsSol(prec,rhs)
%%   function y = subSol(prec,rhs)
    H = prec.H;
    V = prec.V;
    if (size(V,2) == 0)
        y = rhs;
    else
        t = V'*rhs;
        t = H \ t - t;
        y = rhs + V*t;
    end