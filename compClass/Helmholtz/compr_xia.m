function [ A, R, nflops, rk ] = compr_xia( A, tol )

nflops      = 0;
[m,n]       = size(A);
mn          = min(m,n);
piv(1:mn)   = 0;
vn(1:n)     = 0;
R(1:mn,1:n) = 0;
for i = 1:n
    vn(i) = norm(A(1:m,i))^2;
end
nflops = nflops+2*m*n;

rk = mn;
iii = mn;
for i = 1:mn
    [vmax,imax] = max(abs(vn(i:n))); 
    imax = imax+i-1; % important
    if (imax > i)
        A(:,[i imax])     = A(:,[imax i]);
        R(1:i-1,[i imax]) = R(1:i-1,[imax i]);
        piv(i)   = imax;
        vn(imax) = vn(i);
    end
    R(i,i) = norm(A(1:m,i));
    nflops = nflops+2*m;
    
    if (R(i,i) < 1e-15)
        rk  = max(i-1,1);
        iii = i;
        break;
    end

    A(:,i) = A(:,i)/R(i,i);
    nflops = nflops+m;
    vn(i)  = R(i,i);
    R(i,i+1:n) = A(1:m,i)'*A(1:m,i+1:n);
    A(1:m,i+1:n) = A(1:m,i+1:n) - A(1:m,i)*R(i,i+1:n);
    vn(i+1:n) = vn(i+1:n)-abs(R(i,i+1:n)).^2;
    nflops   = nflops+4*m*(n-i);

    if R(i,i) < tol * R(1,1)
        rk  = i;
        iii = i;
        break;
    end
end

A = A(:,1:rk); R = R(1:rk,:);
for i = iii:-1:1
    imax = piv(i);
    mn   = min(rk,rk*imax);
    if mn > 0 
        R(1:mn,[i imax]) = R(1:mn,[imax i]);
    end
end


return;
end