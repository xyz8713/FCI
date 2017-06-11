function  A = E3d_VTI( freq, dim, tos, FS, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, Vp, Vs, rho, epsilon, delta, gamma )
% 3D acoustic Helmholtz equation

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
            
            col_tmp = zeros(3,3,3,3);
            C11 = zeros(3,3,3);
            C12 = zeros(3,3,3);
            C13 = zeros(3,3,3);
            C22 = zeros(3,3,3);
            C23 = zeros(3,3,3);
            C33 = zeros(3,3,3);
            C44 = zeros(3,3,3);
            C55 = zeros(3,3,3);
            C66 = zeros(3,3,3);
            
            for kk = -1:1
                for jj = -1:1
                    for ii = -1:1
                        if xx+ii >= 1 && xx+ii <= Nx && yy+jj >= 1 && yy+jj <= Ny && zz+kk >= 1 && zz+kk <= Nz
                            col_tmp(ii+2,jj+2,kk+2,1) = label(xx+ii,yy+jj,zz+kk)*dim - 2;
                            col_tmp(ii+2,jj+2,kk+2,2) = label(xx+ii,yy+jj,zz+kk)*dim - 1;
                            col_tmp(ii+2,jj+2,kk+2,3) = label(xx+ii,yy+jj,zz+kk)*dim;
                            rr = rho(xx+ii,yy+jj,zz+kk);
                            pp = Vp(xx+ii,yy+jj,zz+kk)^2;
                            ss = Vs(xx+ii,yy+jj,zz+kk)^2;
                            C33(ii+2,jj+2,kk+2) = rr * pp;
                            C11(ii+2,jj+2,kk+2) = C33(ii+2,jj+2,kk+2) * (1 + 2*epsilon(xx+ii,yy+jj,zz+kk));
                            C22(ii+2,jj+2,kk+2) = C11(ii+2,jj+2,kk+2);
                            C44(ii+2,jj+2,kk+2) = rr * ss;
                            C55(ii+2,jj+2,kk+2) = C44(ii+2,jj+2,kk+2);
                            C66(ii+2,jj+2,kk+2) = C44(ii+2,jj+2,kk+2) * (1 + 2*gamma(xx+ii,yy+jj,zz+kk));
                            C12(ii+2,jj+2,kk+2) = C11(ii+2,jj+2,kk+2) - 2*C66(ii+2,jj+2,kk+2);
                            C13(ii+2,jj+2,kk+2) = (C33(ii+2,jj+2,kk+2) - C44(ii+2,jj+2,kk+2)) * sqrt(1+2*delta(xx+ii,yy+jj,zz+kk)/(1-ss/pp)) - C44(ii+2,jj+2,kk+2);
                            C23(ii+2,jj+2,kk+2) = C13(ii+2,jj+2,kk+2);
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%  wavefield one  %%%%%%%%%%%%%%%%%%%%%%%%
            row_tmp = label(xx,yy,zz)*dim - 2;
            val_tmp = zeros(3,3,3,3);
            
            if FS ~= 0 && zz == 1 && C66(2,2,2) ~= 0
                
                val_tmp(2,2,2,1) = -1;
                val_tmp(2,2,3,1) =  1;
                val_tmp(1,2,2,3) = -rx/2 / Sx(xx+tos);
                val_tmp(3,2,2,3) =  rx/2 / Sx(xx+tos);
                
            else
            
                % d^2/dx^2
                tmp  = -rx^2/(C33(2,2,2)*Sx(xx+tos));
                tmp1 = (C11(1,2,2)/Sx(xx-1+tos) + C11(2,2,2)/Sx(xx+tos))/2;
                tmp2 = (C11(3,2,2)/Sx(xx+1+tos) + C11(2,2,2)/Sx(xx+tos))/2;
                val_tmp(1,2,2,1) = val_tmp(1,2,2,1) + tmp * tmp1;
                val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
                val_tmp(3,2,2,1) = val_tmp(3,2,2,1) + tmp * tmp2;

                % d^2/dy^2
                tmp  = -ry^2/(C33(2,2,2)*Sy(yy+tos));
                tmp1 = (C44(2,1,2)/Sy(yy-1+tos) + C44(2,2,2)/Sy(yy+tos))/2;
                tmp2 = (C44(2,3,2)/Sy(yy+1+tos) + C44(2,2,2)/Sy(yy+tos))/2;
                val_tmp(2,1,2,1) = val_tmp(2,1,2,1) + tmp * tmp1;
                val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
                val_tmp(2,3,2,1) = val_tmp(2,3,2,1) + tmp * tmp2;

                % d^2/dz^2
                tmp  = -1/(C33(2,2,2)*Sz(zz+tos));
                tmp1 = (C66(2,2,1)/Sz(zz-1+tos) + C66(2,2,2)/Sz(zz+tos))/2;
                tmp2 = (C66(2,2,3)/Sz(zz+1+tos) + C66(2,2,2)/Sz(zz+tos))/2;
                val_tmp(2,2,1,1) = val_tmp(2,2,1,1) + tmp * tmp1;
                val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - tmp *(tmp1 + tmp2);
                val_tmp(2,2,3,1) = val_tmp(2,2,3,1) + tmp * tmp2;

                % diagonal term
                val_tmp(2,2,2,1) = val_tmp(2,2,2,1) - (hz*freq/Vp(xx,yy,zz))^2;

                % d^2/dxdy
                tmp = -rx*ry/(4*C33(2,2,2)*Sx(xx+tos)*Sy(yy+tos));
                val_tmp(1,1,2,2) = val_tmp(1,1,2,2) + tmp * (C44(2,1,2) + C12(1,2,2));
                val_tmp(3,1,2,2) = val_tmp(3,1,2,2) - tmp * (C44(2,1,2) + C12(3,2,2));
                val_tmp(1,3,2,2) = val_tmp(1,3,2,2) - tmp * (C44(2,3,2) + C12(1,2,2));
                val_tmp(3,3,2,2) = val_tmp(3,3,2,2) + tmp * (C44(2,3,2) + C12(3,2,2));

                % d^2/dxdz
                tmp = -rx/(4*C33(2,2,2)*Sx(xx+tos)*Sz(zz+tos));
                val_tmp(1,2,1,3) = val_tmp(1,2,1,3) + tmp * (C66(2,2,1) + C13(1,2,2));
                val_tmp(3,2,1,3) = val_tmp(3,2,1,3) - tmp * (C66(2,2,1) + C13(3,2,2));
                val_tmp(1,2,3,3) = val_tmp(1,2,3,3) - tmp * (C66(2,2,3) + C13(1,2,2));
                val_tmp(3,2,3,3) = val_tmp(3,2,3,3) + tmp * (C66(2,2,3) + C13(3,2,2));
                
            end

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
            row_tmp = label(xx,yy,zz)*dim - 1;
            val_tmp = zeros(3,3,3,3);
            
            if FS ~= 0 && zz == 1 && C55(2,2,2) ~= 0
                
                val_tmp(2,2,2,2) = -1;
                val_tmp(2,2,3,2) =  1;
                val_tmp(2,1,2,3) = -ry/2 / Sy(yy+tos);
                val_tmp(2,3,2,3) =  ry/2 / Sy(yy+tos);
                
            else
            
                % d^2/dxdy
                tmp = -rx*ry/(4*C33(2,2,2)*Sx(xx+tos)*Sy(yy+tos));
                val_tmp(1,1,2,1) = val_tmp(1,1,2,1) + tmp * (C12(2,1,2) + C44(1,2,2));
                val_tmp(3,1,2,1) = val_tmp(3,1,2,1) - tmp * (C12(2,1,2) + C44(3,2,2));
                val_tmp(1,3,2,1) = val_tmp(1,3,2,1) - tmp * (C12(2,3,2) + C44(1,2,2));
                val_tmp(3,3,2,1) = val_tmp(3,3,2,1) + tmp * (C12(2,3,2) + C44(3,2,2));

                % d^2/dx^2
                tmp  = -rx^2/(C33(2,2,2)*Sx(xx+tos));
                tmp1 = (C44(1,2,2)/Sx(xx-1+tos) + C44(2,2,2)/Sx(xx+tos))/2;
                tmp2 = (C44(3,2,2)/Sx(xx+1+tos) + C44(2,2,2)/Sx(xx+tos))/2;
                val_tmp(1,2,2,2) = val_tmp(1,2,2,2) + tmp * tmp1;
                val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
                val_tmp(3,2,2,2) = val_tmp(3,2,2,2) + tmp * tmp2;

                % d^2/dy^2
                tmp  = -ry^2/(C33(2,2,2)*Sy(yy+tos));
                tmp1 = (C22(2,1,2)/Sy(yy-1+tos) + C22(2,2,2)/Sy(yy+tos))/2;
                tmp2 = (C22(2,3,2)/Sy(yy+1+tos) + C22(2,2,2)/Sy(yy+tos))/2;
                val_tmp(2,1,2,2) = val_tmp(2,1,2,2) + tmp * tmp1;
                val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
                val_tmp(2,3,2,2) = val_tmp(2,3,2,2) + tmp * tmp2;

                % d^2/dz^2
                tmp  = -1/(C33(2,2,2)*Sz(zz+tos));
                tmp1 = (C55(2,2,1)/Sz(zz-1+tos) + C55(2,2,2)/Sz(zz+tos))/2;
                tmp2 = (C55(2,2,3)/Sz(zz+1+tos) + C55(2,2,2)/Sz(zz+tos))/2;
                val_tmp(2,2,1,2) = val_tmp(2,2,1,2) + tmp * tmp1;
                val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - tmp *(tmp1 + tmp2);
                val_tmp(2,2,3,2) = val_tmp(2,2,3,2) + tmp * tmp2;

                % diagonal term
                val_tmp(2,2,2,2) = val_tmp(2,2,2,2) - (hz*freq/Vp(xx,yy,zz))^2;

                % d^2/dydz
                tmp = -ry/(4*C33(2,2,2)*Sy(yy+tos)*Sz(zz+tos));
                val_tmp(2,1,1,3) = val_tmp(2,1,1,3) + tmp * (C55(2,2,1) + C23(2,1,2));
                val_tmp(2,3,1,3) = val_tmp(2,3,1,3) - tmp * (C55(2,2,1) + C23(2,3,2));
                val_tmp(2,1,3,3) = val_tmp(2,1,3,3) - tmp * (C55(2,2,3) + C23(2,1,2));
                val_tmp(2,3,3,3) = val_tmp(2,3,3,3) + tmp * (C55(2,2,3) + C23(2,3,2));
            
            end

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
            
            %%%%%%%%%%%%%%%%%%%%%%%%  wavefield three  %%%%%%%%%%%%%%%%%%%%%%%%
            row_tmp = label(xx,yy,zz)*dim;
            val_tmp = zeros(3,3,3,3);
            
            if FS ~= 0 && zz == 1
                
                val_tmp(1,2,2,1) = -C13(2,2,2)/C33(2,2,2) * rx/2 / Sx(xx+tos);
                val_tmp(3,2,2,1) =  C13(2,2,2)/C33(2,2,2) * rx/2 / Sx(xx+tos);
                val_tmp(2,1,2,2) = -C23(2,2,2)/C33(2,2,2) * ry/2 / Sy(yy+tos);
                val_tmp(2,3,2,2) =  C23(2,2,2)/C33(2,2,2) * ry/2 / Sy(yy+tos);
                val_tmp(2,2,2,3) = -1;
                val_tmp(2,2,3,3) =  1;
                
            else
            
                % d^2/dxdz
                tmp = -rx/(4*C33(2,2,2)*Sx(xx+tos)*Sz(zz+tos));
                val_tmp(1,2,1,1) = val_tmp(1,2,1,1) + tmp * (C13(2,2,1) + C66(1,2,2));
                val_tmp(3,2,1,1) = val_tmp(3,2,1,1) - tmp * (C13(2,2,1) + C66(3,2,2));
                val_tmp(1,2,3,1) = val_tmp(1,2,3,1) - tmp * (C13(2,2,3) + C66(1,2,2));
                val_tmp(3,2,3,1) = val_tmp(3,2,3,1) + tmp * (C13(2,2,3) + C66(3,2,2));

                % d^2/dydz
                tmp = -ry/(4*C33(2,2,2)*Sy(yy+tos)*Sz(zz+tos));
                val_tmp(2,1,1,2) = val_tmp(2,1,1,2) + tmp * (C23(2,2,1) + C55(2,1,2));
                val_tmp(2,3,1,2) = val_tmp(2,3,1,2) - tmp * (C23(2,2,1) + C55(2,3,2));
                val_tmp(2,1,3,2) = val_tmp(2,1,3,2) - tmp * (C23(2,2,3) + C55(2,1,2));
                val_tmp(2,3,3,2) = val_tmp(2,3,3,2) + tmp * (C23(2,2,3) + C55(2,3,2));

                % d^2/dx^2
                tmp  = -rx^2/(C33(2,2,2)*Sx(xx+tos));
                tmp1 = (C66(1,2,2)/Sx(xx-1+tos) + C66(2,2,2)/Sx(xx+tos))/2;
                tmp2 = (C66(3,2,2)/Sx(xx+1+tos) + C66(2,2,2)/Sx(xx+tos))/2;
                val_tmp(1,2,2,3) = val_tmp(1,2,2,3) + tmp * tmp1;
                val_tmp(2,2,2,3) = val_tmp(2,2,2,3) - tmp *(tmp1 + tmp2);
                val_tmp(3,2,2,3) = val_tmp(3,2,2,3) + tmp * tmp2;

                % d^2/dy^2
                tmp  = -ry^2/(C33(2,2,2)*Sy(yy+tos));
                tmp1 = (C55(2,1,2)/Sy(yy-1+tos) + C55(2,2,2)/Sy(yy+tos))/2;
                tmp2 = (C55(2,3,2)/Sy(yy+1+tos) + C55(2,2,2)/Sy(yy+tos))/2;
                val_tmp(2,1,2,3) = val_tmp(2,1,2,3) + tmp * tmp1;
                val_tmp(2,2,2,3) = val_tmp(2,2,2,3) - tmp *(tmp1 + tmp2);
                val_tmp(2,3,2,3) = val_tmp(2,3,2,3) + tmp * tmp2;

                % d^2/dz^2
                tmp  = -1/(C33(2,2,2)*Sz(zz+tos));
                tmp1 = (C33(2,2,1)/Sz(zz-1+tos) + C33(2,2,2)/Sz(zz+tos))/2;
                tmp2 = (C33(2,2,3)/Sz(zz+1+tos) + C33(2,2,2)/Sz(zz+tos))/2;
                val_tmp(2,2,1,3) = val_tmp(2,2,1,3) + tmp * tmp1;
                val_tmp(2,2,2,3) = val_tmp(2,2,2,3) - tmp *(tmp1 + tmp2);
                val_tmp(2,2,3,3) = val_tmp(2,2,3,3) + tmp * tmp2;

                % diagonal term
                val_tmp(2,2,2,3) = val_tmp(2,2,2,3) - (hz*freq/Vp(xx,yy,zz))^2;
            
            end

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