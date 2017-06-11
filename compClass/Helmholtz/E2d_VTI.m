function  A = E2d_VTI( freq, dim, tos, FS, Nz, Nx, N_delta, sigma0, hz, hx, label, Vp, Vs, rho, epsilon, delta, gamma )
% 2D acoustic Helmholtz matrix

rx = hz/hx;

% build up two PML decaying function
Sz = ones(1,Nz+2*tos);
Sx = ones(1,Nx+2*tos);
for ii = 1:N_delta+tos
    tmp = sigma0/sqrt(-1)/freq*(1+cos(pi*ii/(N_delta+tos)))/2;
    Sx(Nx+2*tos-ii+1) = Sx(Nx+2*tos-ii+1) + tmp;
    Sz(Nz+2*tos-ii+1) = Sz(Nz+2*tos-ii+1) + tmp;
    Sx(ii) = Sx(Nx+2*tos-ii+1);
    if FS == 0
        Sz(ii) = Sz(Nz+2*tos-ii+1);
    end 
end

% A generation
len = (9*Nz*Nx - 6*((Nx-2)+(Nz-2)) - 20) * dim^2;
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% No use to follow the tree sequence
for xx = 1:Nx
    for zz = 1:Nz
        
        col_tmp = zeros(3,3,2);
        C11 = zeros(3,3);
        C13 = zeros(3,3);
        C33 = zeros(3,3);
        C66 = zeros(3,3);

        for jj = -1:1
            for ii = -1:1
                if zz+ii >= 1 && zz+ii <= Nz && xx+jj >= 1 && xx+jj <= Nx
                    col_tmp(ii+2,jj+2,1) = label(zz+ii,xx+jj)*dim - 1;
                    col_tmp(ii+2,jj+2,2) = label(zz+ii,xx+jj)*dim;
                    rr = rho(zz+ii,xx+jj);
                    pp = Vp(zz+ii,xx+jj)^2;
                    ss = Vs(zz+ii,xx+jj)^2;
                    C11(ii+2,jj+2) = rr * pp * (1+2*epsilon(zz+ii,xx+jj));
                    C66(ii+2,jj+2) = rr * ss * (1+2*gamma(zz+ii,xx+jj));
                    C33(ii+2,jj+2) = rr * pp;
                    C13(ii+2,jj+2) = rr * ( (pp-ss)*sqrt(1+2*delta(zz+ii,xx+jj)/(1-ss/pp)) - ss );
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%  wavefield one  %%%%%%%%%%%%%%%%%%%%%%%%
        row_tmp = label(zz,xx)*dim - 1;
        val_tmp = zeros(3,3,2);
        
        if FS ~= 0 && zz == 1 && C66(2,2) ~= 0
            
            val_tmp(2,2,1) = -1;
            val_tmp(3,2,1) =  1;
            val_tmp(2,1,2) = -rx/2 / Sx(xx+tos);
            val_tmp(2,3,2) =  rx/2 / Sx(xx+tos);
            
        else
        
            % d^2/dx^2
            tmp  = -rx^2/(C33(2,2)*Sx(xx+tos));
            tmp1 = (C11(2,1)/Sx(xx-1+tos) + C11(2,2)/Sx(xx+tos))/2;
            tmp2 = (C11(2,3)/Sx(xx+1+tos) + C11(2,2)/Sx(xx+tos))/2;
            val_tmp(2,1,1) = val_tmp(2,1,1) + tmp * tmp1;
            val_tmp(2,2,1) = val_tmp(2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,1) = val_tmp(2,3,1) + tmp * tmp2;

            % d^2/dz^2
            tmp = -1/(C33(2,2)*Sz(zz+tos));
            tmp1 = (C66(1,2)/Sz(zz-1+tos) + C66(2,2)/Sz(zz+tos))/2;
            tmp2 = (C66(3,2)/Sz(zz+1+tos) + C66(2,2)/Sz(zz+tos))/2;
            val_tmp(1,2,1) = val_tmp(1,2,1) + tmp * tmp1;
            val_tmp(2,2,1) = val_tmp(2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,1) = val_tmp(3,2,1) + tmp * tmp2;

            % diagonal term
            val_tmp(2,2,1) = val_tmp(2,2,1) - (hz*freq/Vp(zz,xx))^2;

            % d^2/dzdx
            tmp = -rx/(4*C33(2,2)*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,1,2) = val_tmp(1,1,2) + tmp * (C13(2,1) + C66(1,2));
            val_tmp(3,1,2) = val_tmp(3,1,2) - tmp * (C13(2,1) + C66(3,2));
            val_tmp(1,3,2) = val_tmp(1,3,2) - tmp * (C13(2,3) + C66(1,2));
            val_tmp(3,3,2) = val_tmp(3,3,2) + tmp * (C13(2,3) + C66(3,2));
            
        end 

        % storing the values
        for kk = 1:dim
            for jj = 1:3
                for ii = 1:3
                    if col_tmp(ii,jj,kk) ~= 0
                        counter = counter + 1;
                        row(counter) = row_tmp;
                        col(counter) = col_tmp(ii,jj,kk);
                        val(counter) = val_tmp(ii,jj,kk);
                    end
                end
            end 
        end

        %%%%%%%%%%%%%%%%%%%%%%%%  wavefield two  %%%%%%%%%%%%%%%%%%%%%%%%
        row_tmp = label(zz,xx)*dim;
        val_tmp = zeros(3,3,2);
        
        if FS ~= 0 && zz == 1
            
            val_tmp(2,1,1) = -C13(2,2)/C33(2,2) * rx/2 / Sx(xx+tos);
            val_tmp(2,3,1) =  C13(2,2)/C33(2,2) * rx/2 / Sx(xx+tos);
            val_tmp(2,2,2) = -1;
            val_tmp(3,2,2) =  1;
            
        else
        
            % d^2/dzdx
            tmp = -rx/(4*C33(2,2)*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,1,1) = val_tmp(1,1,1) + tmp * (C66(2,1) + C13(1,2));
            val_tmp(3,1,1) = val_tmp(3,1,1) - tmp * (C66(2,1) + C13(3,2));
            val_tmp(1,3,1) = val_tmp(1,3,1) - tmp * (C66(2,3) + C13(1,2));
            val_tmp(3,3,1) = val_tmp(3,3,1) + tmp * (C66(2,3) + C13(3,2));

            % d^2/dx^2
            tmp  = -rx^2/(C33(2,2)*Sx(xx+tos));
            tmp1 = (C66(2,1)/Sx(xx-1+tos) + C66(2,2)/Sx(xx+tos))/2;
            tmp2 = (C66(2,3)/Sx(xx+1+tos) + C66(2,2)/Sx(xx+tos))/2;
            val_tmp(2,1,2) = val_tmp(2,1,2) + tmp * tmp1;
            val_tmp(2,2,2) = val_tmp(2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2) = val_tmp(2,3,2) + tmp * tmp2;

            % d^2/dz^2
            tmp = -1/(C33(2,2)*Sz(zz+tos));
            tmp1 = (C33(1,2)/Sz(zz-1+tos) + C33(2,2)/Sz(zz+tos))/2;
            tmp2 = (C33(3,2)/Sz(zz+1+tos) + C33(2,2)/Sz(zz+tos))/2;
            val_tmp(1,2,2) = val_tmp(1,2,2) + tmp * tmp1;
            val_tmp(2,2,2) = val_tmp(2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2) = val_tmp(3,2,2) + tmp * tmp2;

            % diagonal term
            val_tmp(2,2,2) = val_tmp(2,2,2) - (hz*freq/Vp(zz,xx))^2;
        
        end

        % storing the values
        for kk = 1:dim
            for jj = 1:3
                for ii = 1:3
                    if col_tmp(ii,jj,kk) ~= 0
                        counter = counter + 1;
                        row(counter) = row_tmp;
                        col(counter) = col_tmp(ii,jj,kk);
                        val(counter) = val_tmp(ii,jj,kk);
                    end
                end
            end 
        end
        
    end
end

A = sparse(row,col,val);

return;
end