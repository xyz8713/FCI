     function y = subsSol(prec,rhs)
%%   function y = subSol(prec,rhs)
     H = prec.H;
     V = prec.V;
     t = V'*rhs;
     t = H \ t - t;
     y = rhs + V*t;
