function  A = A3d_CTI( freq, dim, tos, FS, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, Vp, rho, epsilon, delta, theta, phi )
% 3D coupled acoustic TI

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
len = (27*Nx*Ny*Nz - 18*((Nx-2)*(Ny-2)+(Ny-2)*(Nz-2)+(Nz-2)*(Nx-2)) - 60*((Nx-2)+(Ny-2)+(Nz-2)) - 19*8) * dim^2;
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% No use to follow the tree sequence
for zz = 1:Nz
    for yy = 1:Ny
        for xx = 1:Nx

            col_tmp = zeros(3,3,3,2);
            rho_tmp = ones(3,3,3);
            for kk = -1:1
                for jj = -1:1
                    for ii = -1:1
                        if xx+ii >= 1 && xx+ii <= Nx && yy+jj >= 1 && yy+jj <= Ny && zz+kk >= 1 && zz+kk <= Nz
                            col_tmp(ii+2,jj+2,kk+2,1) = label(xx+ii,yy+jj,zz+kk)*dim - 1;
                            col_tmp(ii+2,jj+2,kk+2,2) = label(xx+ii,yy+jj,zz+kk)*dim;
                            rho_tmp(ii+2,jj+2,kk+2)   =   rho(xx+ii,yy+jj,zz+kk);
                        end
                    end
                end
            end
            
            % model parameters
            vv =      Vp(xx,yy,zz);
            ee = epsilon(xx,yy,zz);
            dd =   delta(xx,yy,zz);
            tt =   theta(xx,yy,zz) * pi/180;
            pp =     phi(xx,yy,zz) * pi/180;
            D11 = sin(tt)^2 * cos(pp)^2;
            D12 = sin(tt)^2 * sin(pp)^2;
            D13 = cos(tt)^2;
            D14 = sin(tt)^2 * sin(2*pp);
            D15 = sin(2*tt) * sin(pp);
            D16 = sin(2*tt) * cos(pp);
            D21 = 1 - D11;
            D22 = 1 - D12;
            D23 = 1 - D13;
            D24 = -D14;
            D25 = -D15;
            D26 = -D16;

            %%%%%%%%%%%%%%%%%%%%%%%%  wavefield one  %%%%%%%%%%%%%%%%%%%%%%%%
            row_tmp = label(xx,yy,zz)*dim - 1;
            val_tmp = zeros(3,3,3,2);
            
            %% sigma_H: horizontal stress
            tmp = -(1 + 2*ee);
            Cxx = tmp * D21 * rx^2;
            Cyy = tmp * D22 * ry^2;
            Czz = tmp * D23;
            Cxy = tmp * D24 * rx*ry;
            Cyz = tmp * D25 * ry;
            Czx = tmp * D26 * rx;
            
            % 1. d^2/dx^2
            tmp  = Cxx * rho_tmp(2,2,2) / Sx(xx+tos);
            tmp1 = (1/(rho_tmp(1,2,2)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            tmp2 = (1/(rho_tmp(3,2,2)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            val_tmp(1,2,2,1) = val_tmp(1,2,2,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2,1) = val_tmp(3,2,2,1) + tmp * tmp2;

            % 2. d^2/dy^2
            tmp  = Cyy * rho_tmp(2,2,2) / Sy(yy+tos);
            tmp1 = (1/(rho_tmp(2,1,2)*Sy(yy-1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,3,2)*Sy(yy+1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            val_tmp(2,1,2,1) = val_tmp(2,1,2,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2,1) = val_tmp(2,3,2,1) + tmp * tmp2;

            % 3. d^2/dz^2
            tmp  = Czz * rho_tmp(2,2,2) / Sz(zz+tos);
            tmp1 = (1/(rho_tmp(2,2,1)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,2,3)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            val_tmp(2,2,1,1) = val_tmp(2,2,1,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(2,2,3,1) = val_tmp(2,2,3,1) + tmp * tmp2;

            % 4. d^2/dxdy
            tmp = Cxy * rho_tmp(2,2,2) / (4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(1,1,2,1) = val_tmp(1,1,2,1) + tmp * (rho_tmp(2,1,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(1,2,2));
            val_tmp(3,1,2,1) = val_tmp(3,1,2,1) - tmp * (rho_tmp(2,1,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(3,2,2));
            val_tmp(1,3,2,1) = val_tmp(1,3,2,1) - tmp * (rho_tmp(2,3,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(1,2,2));
            val_tmp(3,3,2,1) = val_tmp(3,3,2,1) + tmp * (rho_tmp(2,3,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(3,2,2));

            % 5. d^2/dydz
            tmp = Cyz * rho_tmp(2,2,2) / (4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(2,1,1,1) = val_tmp(2,1,1,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,1,2));
            val_tmp(2,3,1,1) = val_tmp(2,3,1,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,3,2));
            val_tmp(2,1,3,1) = val_tmp(2,1,3,1) - tmp * (rho_tmp(2,2,3)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,1,2));
            val_tmp(2,3,3,1) = val_tmp(2,3,3,1) + tmp * (rho_tmp(2,2,3)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,3,2));

            % 6. d^2/dzdx
            tmp = Czx * rho_tmp(2,2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,2,1,1) = val_tmp(1,2,1,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(1,2,2));
            val_tmp(3,2,1,1) = val_tmp(3,2,1,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(3,2,2));
            val_tmp(1,2,3,1) = val_tmp(1,2,3,1) - tmp * (rho_tmp(2,2,3)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(1,2,2));
            val_tmp(3,2,3,1) = val_tmp(3,2,3,1) + tmp * (rho_tmp(2,2,3)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(3,2,2));

            % diagonal term
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - (hz*freq/vv)^2;
            
            %% sigma_V: vertical stress
            tmp = -sqrt(1 + 2*dd);
            Cxx = tmp * D11 * rx^2;
            Cyy = tmp * D12 * ry^2;
            Czz = tmp * D13;
            Cxy = tmp * D14 * rx*ry;
            Cyz = tmp * D15 * ry;
            Czx = tmp * D16 * rx;
            
            % 1. d^2/dx^2
            tmp  = Cxx * rho_tmp(2,2,2) / Sx(xx+tos);
            tmp1 = (1/(rho_tmp(1,2,2)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            tmp2 = (1/(rho_tmp(3,2,2)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            val_tmp(1,2,2,2) = val_tmp(1,2,2,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2,2) = val_tmp(3,2,2,2) + tmp * tmp2;

            % 2. d^2/dy^2
            tmp  = Cyy * rho_tmp(2,2,2) / Sy(yy+tos);
            tmp1 = (1/(rho_tmp(2,1,2)*Sy(yy-1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,3,2)*Sy(yy+1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            val_tmp(2,1,2,2) = val_tmp(2,1,2,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2,2) = val_tmp(2,3,2,2) + tmp * tmp2;

            % 3. d^2/dz^2
            tmp  = Czz * rho_tmp(2,2,2) / Sz(zz+tos);
            tmp1 = (1/(rho_tmp(2,2,1)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,2,3)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            val_tmp(2,2,1,2) = val_tmp(2,2,1,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,2,3,2) = val_tmp(2,2,3,2) + tmp * tmp2;

            % 4. d^2/dxdy
            tmp = Cxy * rho_tmp(2,2,2) / (4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(1,1,2,2) = val_tmp(1,1,2,2) + tmp * (rho_tmp(2,1,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(1,2,2));
            val_tmp(3,1,2,2) = val_tmp(3,1,2,2) - tmp * (rho_tmp(2,1,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(3,2,2));
            val_tmp(1,3,2,2) = val_tmp(1,3,2,2) - tmp * (rho_tmp(2,3,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(1,2,2));
            val_tmp(3,3,2,2) = val_tmp(3,3,2,2) + tmp * (rho_tmp(2,3,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(3,2,2));

            % 5. d^2/dydz
            tmp = Cyz * rho_tmp(2,2,2) / (4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(2,1,1,2) = val_tmp(2,1,1,2) + tmp * (rho_tmp(2,2,1)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,1,2));
            val_tmp(2,3,1,2) = val_tmp(2,3,1,2) - tmp * (rho_tmp(2,2,1)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,3,2));
            val_tmp(2,1,3,2) = val_tmp(2,1,3,2) - tmp * (rho_tmp(2,2,3)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,1,2));
            val_tmp(2,3,3,2) = val_tmp(2,3,3,2) + tmp * (rho_tmp(2,2,3)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,3,2));

            % 6. d^2/dzdx
            tmp = Czx * rho_tmp(2,2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,2,1,2) = val_tmp(1,2,1,2) + tmp * (rho_tmp(2,2,1)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(1,2,2));
            val_tmp(3,2,1,2) = val_tmp(3,2,1,2) - tmp * (rho_tmp(2,2,1)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(3,2,2));
            val_tmp(1,2,3,2) = val_tmp(1,2,3,2) - tmp * (rho_tmp(2,2,3)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(1,2,2));
            val_tmp(3,2,3,2) = val_tmp(3,2,3,2) + tmp * (rho_tmp(2,2,3)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(3,2,2));

            % storing the values
            for mm = 1:dim
                for kk = 1:3
                    for jj = 1:3
                        for ii = 1:3
                            if col_tmp(ii,jj,kk,mm) ~= 0
                                counter = counter + 1;
                                row(counter) = row_tmp;
                                col(counter) = col_tmp(ii,jj,kk,mm);
                                val(counter) = val_tmp(ii,jj,kk,mm);
                            end
                        end
                    end
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%  wavefield two  %%%%%%%%%%%%%%%%%%%%%%%%
            row_tmp = label(xx,yy,zz)*dim;
            val_tmp = zeros(3,3,3,2);
            
            %% sigma_H: horizontal stress
            tmp = -sqrt(1 + 2.0*dd);
            Cxx = tmp * D21 * rx^2;
            Cyy = tmp * D22 * ry^2;
            Czz = tmp * D23;
            Cxy = tmp * D24 * rx*ry;
            Cyz = tmp * D25 * ry;
            Czx = tmp * D26 * rx;
            
            % 1. d^2/dx^2
            tmp  = Cxx * rho_tmp(2,2,2) / Sx(xx+tos);
            tmp1 = (1/(rho_tmp(1,2,2)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            tmp2 = (1/(rho_tmp(3,2,2)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            val_tmp(1,2,2,1) = val_tmp(1,2,2,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2,1) = val_tmp(3,2,2,1) + tmp * tmp2;

            % 2. d^2/dy^2
            tmp  = Cyy * rho_tmp(2,2,2) / Sy(yy+tos);
            tmp1 = (1/(rho_tmp(2,1,2)*Sy(yy-1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,3,2)*Sy(yy+1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            val_tmp(2,1,2,1) = val_tmp(2,1,2,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2,1) = val_tmp(2,3,2,1) + tmp * tmp2;

            % 3. d^2/dz^2
            tmp  = Czz * rho_tmp(2,2,2) / Sz(zz+tos);
            tmp1 = (1/(rho_tmp(2,2,1)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,2,3)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            val_tmp(2,2,1,1) = val_tmp(2,2,1,1) + tmp * tmp1;
            val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
            val_tmp(2,2,3,1) = val_tmp(2,2,3,1) + tmp * tmp2;

            % 4. d^2/dxdy
            tmp = Cxy * rho_tmp(2,2,2) / (4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(1,1,2,1) = val_tmp(1,1,2,1) + tmp * (rho_tmp(2,1,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(1,2,2));
            val_tmp(3,1,2,1) = val_tmp(3,1,2,1) - tmp * (rho_tmp(2,1,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(3,2,2));
            val_tmp(1,3,2,1) = val_tmp(1,3,2,1) - tmp * (rho_tmp(2,3,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(1,2,2));
            val_tmp(3,3,2,1) = val_tmp(3,3,2,1) + tmp * (rho_tmp(2,3,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(3,2,2));

            % 5. d^2/dydz
            tmp = Cyz * rho_tmp(2,2,2) / (4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(2,1,1,1) = val_tmp(2,1,1,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,1,2));
            val_tmp(2,3,1,1) = val_tmp(2,3,1,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,3,2));
            val_tmp(2,1,3,1) = val_tmp(2,1,3,1) - tmp * (rho_tmp(2,2,3)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,1,2));
            val_tmp(2,3,3,1) = val_tmp(2,3,3,1) + tmp * (rho_tmp(2,2,3)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,3,2));

            % 6. d^2/dzdx
            tmp = Czx * rho_tmp(2,2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,2,1,1) = val_tmp(1,2,1,1) + tmp * (rho_tmp(2,2,1)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(1,2,2));
            val_tmp(3,2,1,1) = val_tmp(3,2,1,1) - tmp * (rho_tmp(2,2,1)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(3,2,2));
            val_tmp(1,2,3,1) = val_tmp(1,2,3,1) - tmp * (rho_tmp(2,2,3)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(1,2,2));
            val_tmp(3,2,3,1) = val_tmp(3,2,3,1) + tmp * (rho_tmp(2,2,3)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(3,2,2));
            
            %% sigma_V: vertical stress
            Cxx = -D11 * rx^2;
            Cyy = -D12 * ry^2;
            Czz = -D13;
            Cxy = -D14 * rx*ry;
            Cyz = -D15 * ry;
            Czx = -D16 * rx;
            
            % 1. d^2/dx^2
            tmp  = Cxx * rho_tmp(2,2,2) / Sx(xx+tos);
            tmp1 = (1/(rho_tmp(1,2,2)*Sx(xx-1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            tmp2 = (1/(rho_tmp(3,2,2)*Sx(xx+1+tos)) + 1/(rho_tmp(2,2,2)*Sx(xx+tos))) / 2;
            val_tmp(1,2,2,2) = val_tmp(1,2,2,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(3,2,2,2) = val_tmp(3,2,2,2) + tmp * tmp2;

            % 2. d^2/dy^2
            tmp  = Cyy * rho_tmp(2,2,2) / Sy(yy+tos);
            tmp1 = (1/(rho_tmp(2,1,2)*Sy(yy-1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,3,2)*Sy(yy+1+tos)) + 1/(rho_tmp(2,2,2)*Sy(yy+tos))) / 2;
            val_tmp(2,1,2,2) = val_tmp(2,1,2,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,3,2,2) = val_tmp(2,3,2,2) + tmp * tmp2;

            % 3. d^2/dz^2
            tmp  = Czz * rho_tmp(2,2,2) / Sz(zz+tos);
            tmp1 = (1/(rho_tmp(2,2,1)*Sz(zz-1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            tmp2 = (1/(rho_tmp(2,2,3)*Sz(zz+1+tos)) + 1/(rho_tmp(2,2,2)*Sz(zz+tos))) / 2;
            val_tmp(2,2,1,2) = val_tmp(2,2,1,2) + tmp * tmp1;
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
            val_tmp(2,2,3,2) = val_tmp(2,2,3,2) + tmp * tmp2;

            % 4. d^2/dxdy
            tmp = Cxy * rho_tmp(2,2,2) / (4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(1,1,2,2) = val_tmp(1,1,2,2) + tmp * (rho_tmp(2,1,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(1,2,2));
            val_tmp(3,1,2,2) = val_tmp(3,1,2,2) - tmp * (rho_tmp(2,1,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,1,2)*rho_tmp(3,2,2));
            val_tmp(1,3,2,2) = val_tmp(1,3,2,2) - tmp * (rho_tmp(2,3,2)+rho_tmp(1,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(1,2,2));
            val_tmp(3,3,2,2) = val_tmp(3,3,2,2) + tmp * (rho_tmp(2,3,2)+rho_tmp(3,2,2))/(2*rho_tmp(2,3,2)*rho_tmp(3,2,2));

            % 5. d^2/dydz
            tmp = Cyz * rho_tmp(2,2,2) / (4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(2,1,1,2) = val_tmp(2,1,1,2) + tmp * (rho_tmp(2,2,1)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,1,2));
            val_tmp(2,3,1,2) = val_tmp(2,3,1,2) - tmp * (rho_tmp(2,2,1)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,1)*rho_tmp(2,3,2));
            val_tmp(2,1,3,2) = val_tmp(2,1,3,2) - tmp * (rho_tmp(2,2,3)+rho_tmp(2,1,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,1,2));
            val_tmp(2,3,3,2) = val_tmp(2,3,3,2) + tmp * (rho_tmp(2,2,3)+rho_tmp(2,3,2))/(2*rho_tmp(2,2,3)*rho_tmp(2,3,2));

            % 6. d^2/dzdx
            tmp = Czx * rho_tmp(2,2,2) / (4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(1,2,1,2) = val_tmp(1,2,1,2) + tmp * (rho_tmp(2,2,1)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(1,2,2));
            val_tmp(3,2,1,2) = val_tmp(3,2,1,2) - tmp * (rho_tmp(2,2,1)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,1)*rho_tmp(3,2,2));
            val_tmp(1,2,3,2) = val_tmp(1,2,3,2) - tmp * (rho_tmp(2,2,3)+rho_tmp(1,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(1,2,2));
            val_tmp(3,2,3,2) = val_tmp(3,2,3,2) + tmp * (rho_tmp(2,2,3)+rho_tmp(3,2,2))/(2*rho_tmp(2,2,3)*rho_tmp(3,2,2));
            
            % 7. diagonal term
            val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - (hz*freq/vv)^2;

            % storing the values
            for mm = 1:dim
                for kk = 1:3
                    for jj = 1:3
                        for ii = 1:3
                            if col_tmp(ii,jj,kk,mm) ~= 0
                                counter = counter + 1;
                                row(counter) = row_tmp;
                                col(counter) = col_tmp(ii,jj,kk,mm);
                                val(counter) = val_tmp(ii,jj,kk,mm);
                            end
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