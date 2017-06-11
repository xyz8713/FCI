   function [D] = ExEigLap3(nx, ny, nz, a, b, shift)
%% function [D] = ExEigLap3(nx, ny, nz, a, b, shift)
%  returns all eigenvalues of 2-D laplacean in [a, b] 
if (nargin <6) 
    shift = 0;
end
a = a+shift;
b = b+shift;
thetax = 0.5*pi/(nx+1);
thetay = 0.5*pi/(ny+1);
thetaz = 0.5*pi/(nz+1);
l = 0 ;
D = []; 
for i=1:nx
    tx = sin(i*thetax) ;
    for j=1:ny
        ty = sin(j*thetay);
        for k=1:nz
          tz = sin(k*thetaz);
          if (nz==1), tz=0;, end
          eval = 4*(tx*tx + ty*ty + tz*tz);
            if (eval >= a & eval <= b)
                l = l+1;
                D(l) = eval-shift;
            end
        end
    end
end
D = sort(D);
D = D(:);