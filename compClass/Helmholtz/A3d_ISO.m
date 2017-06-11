function  A = A3d_ISO( freq, FS, Nx, Ny, Nz, N_delta, sigma0, label, vel)
% 3D acoustic Helmholtz equation

% simplified by Xiao
hx = 1/Nx;
hy = 1/Ny;
hz = 1/Nz;
tos = 1;
rx = hz/hx;
ry = hz/hy;

% build up three PML decaying function
Sx = ones(1,Nx+2*tos);
Sy = ones(1,Ny+2*tos);
Sz = ones(1,Nz+2*tos);
for ii = 1:N_delta+tos
    tmp = sigma0/sqrt(-1)/freq*(1+cos(pi*ii/(N_delta+tos)))/2;
    Sx(Nx+2*tos-ii+1) = Sx(Nx+2*tos-ii+1) + tmp;
    Sy(Ny+2*tos-ii+1) = Sy(Ny+2*tos-ii+1) + tmp;
    Sz(Nz+2*tos-ii+1) = Sz(Nz+2*tos-ii+1) + tmp;
    Sx(ii) = Sx(Nx+2*tos-ii+1);
    Sy(ii) = Sy(Ny+2*tos-ii+1);
    if FS == 0
        Sz(ii) = Sz(Nz+2*tos-ii+1);
    end 
end

% A generation
len = 27*Nx*Ny*Nz - 18*((Nx-2)*(Ny-2)+(Ny-2)*(Nz-2)+(Nz-2)*(Nx-2)) - 60*((Nx-2)+(Ny-2)+(Nz-2)) - 19*8;
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% No use to follow the tree sequence
for zz = 1:Nz
    for yy = 1:Ny
        for xx = 1:Nx

            row_tmp = label(xx,yy,zz) * ones(3,3,3);
            col_tmp = zeros(3,3,3);
            rho_tmp = ones(3,3,3);
            for kk = -1:1
                for jj = -1:1
                    for ii = -1:1
                        if xx+ii >= 1 && xx+ii <= Nx && yy+jj >= 1 && yy+jj <= Ny && zz+kk >= 1 && zz+kk <= Nz
                            col_tmp(ii+2,jj+2,kk+2) = label(xx+ii,yy+jj,zz+kk);
                        end
                    end
                end
            end

            % value preprocessing
            Cxx = -1 * rx^2;
            Cyy = -1 * ry^2;
            Czz = -1;
            Cxy = 0 * rx*ry;
            Cyz = 0 * ry;
            Czx = 0 * rx;
            
            % value generating
            val_tmp = zeros(3,3,3);
            
            % 1. d^2/dx^2
            tmp  = Cxx * rho_tmp(2,2,2) / Sx(xx+tos);
            tmp1 = (1/(rho_tmp(1,2,2)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos)))/2;
            tmp2 = (1/(rho_tmp(3,2,2)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos)))/2;
            val_tmp(1,2,2) = val_tmp(1,2,2) + tmp * tmp1;
            val_tmp(2,2,2) = val_tmp(2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2) = val_tmp(3,2,2) + tmp * tmp2;

            % 2. d^2/dy^2
            tmp  = Cyy * rho_tmp(2,2,2) / Sy(yy+tos);
            tmp1 = (1/(rho_tmp(2,1,2)*Sy(yy-1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos)))/2;
            tmp2 = (1/(rho_tmp(2,3,2)*Sy(yy+1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos)))/2;
            val_tmp(2,1,2) = val_tmp(2,1,2) + tmp * tmp1;
            val_tmp(2,2,2) = val_tmp(2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2) = val_tmp(2,3,2) + tmp * tmp2;

            % 3. d^2/dz^2
            tmp  = Czz * rho_tmp(2,2,2) / Sz(zz+tos);
            tmp1 = (1/(rho_tmp(2,2,1)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos)))/2;
            tmp2 = (1/(rho_tmp(2,2,3)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos)))/2;
            val_tmp(2,2,1) = val_tmp(2,2,1) + tmp * tmp1;
            val_tmp(2,2,2) = val_tmp(2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,2,3) = val_tmp(2,2,3) + tmp * tmp2;

            % 4. d^2/dxdy
            tmp = Cxy * rho_tmp(2,2,2) / (4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(1,1,2) = val_tmp(1,1,2) + tmp * (rho_tmp(2,1,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(1,2,2));
            val_tmp(3,1,2) = val_tmp(3,1,2) - tmp * (rho_tmp(2,1,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(3,2,2));
            val_tmp(1,3,2) = val_tmp(1,3,2) - tmp * (rho_tmp(2,3,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(1,2,2));
            val_tmp(3,3,2) = val_tmp(3,3,2) + tmp * (rho_tmp(2,3,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(3,2,2));

            % 5. d^2/dydz
            tmp = Cyz * rho_tmp(2,2,2) / (4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(2,1,1) = val_tmp(2,1,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,1,2));
            val_tmp(2,3,1) = val_tmp(2,3,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,3,2));
            val_tmp(2,1,3) = val_tmp(2,1,3) - tmp * (rho_tmp(2,2,3)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,1,2));
            val_tmp(2,3,3) = val_tmp(2,3,3) + tmp * (rho_tmp(2,2,3)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,3,2));

            % 6. d^2/dzdx
            tmp = Czx * rho_tmp(2,2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,2,1) = val_tmp(1,2,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(1,2,2));
            val_tmp(3,2,1) = val_tmp(3,2,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(3,2,2));
            val_tmp(1,2,3) = val_tmp(1,2,3) - tmp * (rho_tmp(2,2,3)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(1,2,2));
            val_tmp(3,2,3) = val_tmp(3,2,3) + tmp * (rho_tmp(2,2,3)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(3,2,2));
            
            % mass term
            val_tmp(2,2,2) = val_tmp(2,2,2) - (hz*freq/vel(xx,yy,zz))^2;
            
            for kk = 1:3
                for jj = 1:3
                    for ii = 1:3
                        if col_tmp(ii,jj,kk) ~= 0
                            counter = counter+1;
                            row(counter) = row_tmp(ii,jj,kk);
                            col(counter) = col_tmp(ii,jj,kk);
                            val(counter) = val_tmp(ii,jj,kk);
                        end
                    end
                end
            end

        end
    end
end

A = sparse(row,col,val);

return;
end