     function y = subsSol(prec,rhs)
%%   function y = subSol(prec,rhs)
     H = prec.H;
     V = prec.V;
     U = prec.U;
     t = U'*rhs;
     y = V*(H\t);
