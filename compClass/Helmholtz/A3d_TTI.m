function  A = A3d_TTI( freq, tos, FS, Nx, Ny, Nz, N_delta, sigma0, hx, hy, hz, label, vel_pz, rho, epsilon, delta, theta, phi )

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
len = 125*Nx*Ny*Nz - 8*((125-27)+(125-36)*3+(125-48)*3+(125-64)) - 4*((125-45)+(125-60)*2+(125-80))*((Nx-4)+(Ny-4)+(Nz-4)) ...
      - 2*(25+50)*((Nx-4)*(Ny-4)+(Ny-4)*(Nz-4)+(Nz-4)*(Nx-4));
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% main loops
for zz = 1:Nz
    for yy = 1:Ny
        for xx = 1:Nx

            row_tmp = label(xx,yy,zz) * ones(5,5,5);
            col_tmp = zeros(5,5,5);
            rho_tmp = ones(5,5,5);
            
            for kk = -2:2
                for jj = -2:2
                    for ii = -2:2
                        if xx+ii >= 1 && xx+ii <= Nx && yy+jj >= 1 && yy+jj <= Ny && zz+kk >= 1 && zz+kk <= Nz
                            col_tmp(ii+3,jj+3,kk+3) = label(xx+ii,yy+jj,zz+kk);
                            rho_tmp(ii+3,jj+3,kk+3) =   rho(xx+ii,yy+jj,zz+kk);
                        end
                    end
                end
            end
            
            vp  =  vel_pz(xx,yy,zz);
            ee  = epsilon(xx,yy,zz);
            dd  =   delta(xx,yy,zz);
            tt  =   theta(xx,yy,zz)*pi/180;
            pp  =     phi(xx,yy,zz)*pi/180;
            D11 = sin(tt)^2*cos(pp)^2;    D12 = sin(tt)^2*sin(pp)^2;    D13 = cos(tt)^2;
            D14 = sin(tt)^2*sin(2*pp);    D15 = sin(2*tt)*sin(pp);      D16 = sin(2*tt)*cos(pp);
            D21 = 1 - D11;                D22 = 1 - D12;                D23 = 1 - D13;
            D24 = -D14;                   D25 = -D15;                   D26 = -D16;
            
            Cxx   = (-D11 - (1+2*ee)*D21) * rx^2;
            Cyy   = (-D12 - (1+2*ee)*D22) * ry^2;
            Czz   =  -D13 - (1+2*ee)*D23;
            Cxy   = (-D14 - (1+2*ee)*D24) * rx*ry;
            Cyz   = (-D15 - (1+2*ee)*D25) * ry;
            Czx   = (-D16 - (1+2*ee)*D26) * rx;
            
            Cxxxx = (-2*(ee-dd)*(vp/hz/freq)^2*D11*D21) * rx^4;
            Cyyyy = (-2*(ee-dd)*(vp/hz/freq)^2*D12*D22) * ry^4;
            Czzzz =  -2*(ee-dd)*(vp/hz/freq)^2*D13*D23;
            
            Cxxxy = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D24+D14*D21)) * rx^3*ry;
            Cxyyy = (-2*(ee-dd)*(vp/hz/freq)^2*(D12*D24+D14*D22)) * ry^3*rx;
            Cxxxz = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D26+D16*D21)) * rx^3;
            Cxzzz = (-2*(ee-dd)*(vp/hz/freq)^2*(D13*D26+D16*D23)) * rx;
            Cyyyz = (-2*(ee-dd)*(vp/hz/freq)^2*(D12*D25+D15*D22)) * ry^3;
            Cyzzz = (-2*(ee-dd)*(vp/hz/freq)^2*(D13*D25+D15*D23)) * ry;
            Cxxyz = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D25+D15*D21)) * rx^2*ry;
            Cxyyz = (-2*(ee-dd)*(vp/hz/freq)^2*(D12*D26+D16*D22)) * ry^2*rx;
            Cxyzz = (-2*(ee-dd)*(vp/hz/freq)^2*(D13*D24+D14*D23)) * rx*ry;
            
            Cxxyy = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D22+D12*D21+D14*D24)) * rx^2*ry^2;
            Cyyzz = (-2*(ee-dd)*(vp/hz/freq)^2*(D12*D23+D13*D22+D15*D25)) * ry^2;
            Cxxzz = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D23+D13*D21+D16*D26)) * rx^2;
            
            
            % value generating
            val_tmp = zeros(5,5,5);
            
            
            %%%%%%%%%%%%%%%%%%%%%  1. x^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxx * rho_tmp(3,3,3)/Sx(xx+tos);
            tmp1 = (1/(rho_tmp(2,3,3)*Sx(xx-1+tos)) + 1/(rho_tmp(3,3,3)*Sx(xx+tos)))/2;
            tmp2 = (1/(rho_tmp(4,3,3)*Sx(xx+1+tos)) + 1/(rho_tmp(3,3,3)*Sx(xx+tos)))/2;
            val_tmp(2,3,3) = val_tmp(2,3,3) + tmp * tmp1;
            val_tmp(3,3,3) = val_tmp(3,3,3) - tmp *(tmp1 + tmp2);
            val_tmp(4,3,3) = val_tmp(4,3,3) + tmp * tmp2;
            
            
            %%%%%%%%%%%%%%%%%%%%%  2. y^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyy * rho_tmp(3,3,3)/Sy(yy+tos);
            tmp1 = (1/(rho_tmp(3,2,3)*Sy(yy-1+tos)) + 1/(rho_tmp(3,3,3)*Sy(yy+tos)))/2;
            tmp2 = (1/(rho_tmp(3,4,3)*Sy(yy+1+tos)) + 1/(rho_tmp(3,3,3)*Sy(yy+tos)))/2;
            val_tmp(3,2,3) = val_tmp(3,2,3) + tmp * tmp1;
            val_tmp(3,3,3) = val_tmp(3,3,3) - tmp *(tmp1 + tmp2);
            val_tmp(3,4,3) = val_tmp(3,4,3) + tmp * tmp2;
            
            
            %%%%%%%%%%%%%%%%%%%%%  3. z^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Czz * rho_tmp(3,3,3)/Sz(zz+tos);
            tmp1 = (1/(rho_tmp(3,3,2)*Sz(zz-1+tos)) + 1/(rho_tmp(3,3,3)*Sz(zz+tos)))/2;
            tmp2 = (1/(rho_tmp(3,3,4)*Sz(zz+1+tos)) + 1/(rho_tmp(3,3,3)*Sz(zz+tos)))/2;
            val_tmp(3,3,2) = val_tmp(3,3,2) + tmp * tmp1;
            val_tmp(3,3,3) = val_tmp(3,3,3) - tmp *(tmp1 + tmp2);
            val_tmp(3,3,4) = val_tmp(3,3,4) + tmp * tmp2;
            
            
            %%%%%%%%%%%%%%%%%%%%%  4. xy   %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxy * rho_tmp(3,3,3)/(4*Sx(xx+tos)*Sy(yy+tos));
            val_tmp(2,2,3) = val_tmp(2,2,3) + tmp * (rho_tmp(3,2,3)+rho_tmp(2,3,3))/(2*rho_tmp(3,2,3)*rho_tmp(2,3,3));
            val_tmp(4,2,3) = val_tmp(4,2,3) - tmp * (rho_tmp(3,2,3)+rho_tmp(4,3,3))/(2*rho_tmp(3,2,3)*rho_tmp(4,3,3));
            val_tmp(2,4,3) = val_tmp(2,4,3) - tmp * (rho_tmp(3,4,3)+rho_tmp(2,3,3))/(2*rho_tmp(3,4,3)*rho_tmp(2,3,3));
            val_tmp(4,4,3) = val_tmp(4,4,3) + tmp * (rho_tmp(3,4,3)+rho_tmp(4,3,3))/(2*rho_tmp(3,4,3)*rho_tmp(4,3,3));
            
            
            %%%%%%%%%%%%%%%%%%%%%  5. yz   %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyz * rho_tmp(3,3,3)/(4*Sy(yy+tos)*Sz(zz+tos));
            val_tmp(3,2,2) = val_tmp(3,2,2) + tmp * (rho_tmp(3,3,2)+rho_tmp(3,2,3))/(2*rho_tmp(3,3,2)*rho_tmp(3,2,3));
            val_tmp(3,4,2) = val_tmp(3,4,2) - tmp * (rho_tmp(3,3,2)+rho_tmp(3,4,3))/(2*rho_tmp(3,3,2)*rho_tmp(3,4,3));
            val_tmp(3,2,4) = val_tmp(3,2,4) - tmp * (rho_tmp(3,2,3)+rho_tmp(3,3,4))/(2*rho_tmp(3,2,3)*rho_tmp(3,3,4));
            val_tmp(3,4,4) = val_tmp(3,4,4) + tmp * (rho_tmp(3,4,3)+rho_tmp(3,3,4))/(2*rho_tmp(3,4,3)*rho_tmp(3,3,4));
         
            
            %%%%%%%%%%%%%%%%%%%%%  6. zx   %%%%%%%%%%%%%%%%%%%%%%%%%% 
            tmp = Czx * rho_tmp(3,3,3)/(4*Sz(zz+tos)*Sx(xx+tos));
            val_tmp(2,3,2) = val_tmp(2,3,2) + tmp * (rho_tmp(3,3,2)+rho_tmp(2,3,3))/(2*rho_tmp(3,3,2)*rho_tmp(2,3,3));
            val_tmp(4,3,2) = val_tmp(4,3,2) - tmp * (rho_tmp(3,3,2)+rho_tmp(4,3,3))/(2*rho_tmp(3,3,2)*rho_tmp(4,3,3));
            val_tmp(2,3,4) = val_tmp(2,3,4) - tmp * (rho_tmp(2,3,3)+rho_tmp(3,3,4))/(2*rho_tmp(2,3,3)*rho_tmp(3,3,4));
            val_tmp(4,3,4) = val_tmp(4,3,4) + tmp * (rho_tmp(4,3,3)+rho_tmp(3,3,4))/(2*rho_tmp(4,3,3)*rho_tmp(3,3,4));
            
            
            %%%%%%%%%%%%%%%%%%%%%  7. x^4  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxxx / Sx(xx+tos)^4;
            val_tmp(1,3,3) = val_tmp(1,3,3) + tmp;
            val_tmp(2,3,3) = val_tmp(2,3,3) - tmp*4;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*6;
            val_tmp(4,3,3) = val_tmp(4,3,3) - tmp*4;
            val_tmp(5,3,3) = val_tmp(5,3,3) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%%  8. y^4  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyyyy / Sy(yy+tos)^4;
            val_tmp(3,1,3) = val_tmp(3,1,3) + tmp;
            val_tmp(3,2,3) = val_tmp(3,2,3) - tmp*4;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*6;
            val_tmp(3,4,3) = val_tmp(3,4,3) - tmp*4;
            val_tmp(3,5,3) = val_tmp(3,5,3) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%%  9. z^4  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Czzzz / Sz(zz+tos)^4;
            val_tmp(3,3,1) = val_tmp(3,3,1) + tmp;
            val_tmp(3,3,2) = val_tmp(3,3,2) - tmp*4;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*6;
            val_tmp(3,3,4) = val_tmp(3,3,4) - tmp*4;
            val_tmp(3,3,5) = val_tmp(3,3,5) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%% 10. x^3 y  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxxy / (Sx(xx+tos)^3 * Sy(yy+tos));
            val_tmp(1,2,3) = val_tmp(1,2,3) + tmp /4;
            val_tmp(2,2,3) = val_tmp(2,2,3) - tmp /2;
            val_tmp(4,2,3) = val_tmp(4,2,3) + tmp /2;
            val_tmp(5,2,3) = val_tmp(5,2,3) - tmp /4;
            val_tmp(1,4,3) = val_tmp(1,4,3) - tmp /4;
            val_tmp(2,4,3) = val_tmp(2,4,3) + tmp /2;
            val_tmp(4,4,3) = val_tmp(4,4,3) - tmp /2;
            val_tmp(5,4,3) = val_tmp(5,4,3) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 11. x y^3  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxyyy / (Sx(xx+tos) * Sy(yy+tos)^3);
            val_tmp(2,1,3) = val_tmp(2,1,3) + tmp /4;
            val_tmp(2,2,3) = val_tmp(2,2,3) - tmp /2;
            val_tmp(2,4,3) = val_tmp(2,4,3) + tmp /2;
            val_tmp(2,5,3) = val_tmp(2,5,3) - tmp /4;
            val_tmp(4,1,3) = val_tmp(4,1,3) - tmp /4;
            val_tmp(4,2,3) = val_tmp(4,2,3) + tmp /2;
            val_tmp(4,4,3) = val_tmp(4,4,3) - tmp /2;
            val_tmp(4,5,3) = val_tmp(4,5,3) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 12. x^3 z  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxxz / (Sx(xx+tos)^3 * Sz(zz+tos));
            val_tmp(1,3,2) = val_tmp(1,3,2) + tmp /4;
            val_tmp(2,3,2) = val_tmp(2,3,2) - tmp /2;
            val_tmp(4,3,2) = val_tmp(4,3,2) + tmp /2;
            val_tmp(5,3,2) = val_tmp(5,3,2) - tmp /4;
            val_tmp(1,3,4) = val_tmp(1,3,4) - tmp /4;
            val_tmp(2,3,4) = val_tmp(2,3,4) + tmp /2;
            val_tmp(4,3,4) = val_tmp(4,3,4) - tmp /2;
            val_tmp(5,3,4) = val_tmp(5,3,4) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 13. x z^3  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxzzz / (Sx(xx+tos) * Sz(zz+tos)^3);
            val_tmp(2,3,1) = val_tmp(2,3,1) + tmp /4;
            val_tmp(2,3,2) = val_tmp(2,3,2) - tmp /2;
            val_tmp(2,3,4) = val_tmp(2,3,4) + tmp /2;
            val_tmp(2,3,5) = val_tmp(2,3,5) - tmp /4;
            val_tmp(4,3,1) = val_tmp(4,3,1) - tmp /4;
            val_tmp(4,3,2) = val_tmp(4,3,2) + tmp /2;
            val_tmp(4,3,4) = val_tmp(4,3,4) - tmp /2;
            val_tmp(4,3,5) = val_tmp(4,3,5) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 14. y^3 z  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyyyz / (Sy(yy+tos)^3 * Sz(zz+tos));
            val_tmp(3,1,2) = val_tmp(3,1,2) + tmp /4;
            val_tmp(3,2,2) = val_tmp(3,2,2) - tmp /2;
            val_tmp(3,4,2) = val_tmp(3,4,2) + tmp /2;
            val_tmp(3,5,2) = val_tmp(3,5,2) - tmp /4;
            val_tmp(3,1,4) = val_tmp(3,1,4) - tmp /4;
            val_tmp(3,2,4) = val_tmp(3,2,4) + tmp /2;
            val_tmp(3,4,4) = val_tmp(3,4,4) - tmp /2;
            val_tmp(3,5,4) = val_tmp(3,5,4) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 15. y z^3  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyzzz / (Sy(yy+tos) * Sz(zz+tos)^3);
            val_tmp(3,2,1) = val_tmp(3,2,1) + tmp /4;
            val_tmp(3,2,2) = val_tmp(3,2,2) - tmp /2;
            val_tmp(3,2,4) = val_tmp(3,2,4) + tmp /2;
            val_tmp(3,2,5) = val_tmp(3,2,5) - tmp /4;
            val_tmp(3,4,1) = val_tmp(3,4,1) - tmp /4;
            val_tmp(3,4,2) = val_tmp(3,4,2) + tmp /2;
            val_tmp(3,4,4) = val_tmp(3,4,4) - tmp /2;
            val_tmp(3,4,5) = val_tmp(3,4,5) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 16. x^2 yz  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxyz / (Sx(xx+tos)^2 * Sy(yy+tos) * Sz(zz+tos));
            val_tmp(2,2,2) = val_tmp(2,2,2) + tmp /4;
            val_tmp(2,4,2) = val_tmp(2,4,2) - tmp /4;
            val_tmp(2,2,4) = val_tmp(2,2,4) - tmp /4;
            val_tmp(2,4,4) = val_tmp(2,4,4) + tmp /4;
            val_tmp(3,2,2) = val_tmp(3,2,2) - tmp /2;
            val_tmp(3,4,2) = val_tmp(3,4,2) + tmp /2;
            val_tmp(3,2,4) = val_tmp(3,2,4) + tmp /2;
            val_tmp(3,4,4) = val_tmp(3,4,4) - tmp /2;
            val_tmp(4,2,2) = val_tmp(4,2,2) + tmp /4;
            val_tmp(4,4,2) = val_tmp(4,4,2) - tmp /4;
            val_tmp(4,2,4) = val_tmp(4,2,4) - tmp /4;
            val_tmp(4,4,4) = val_tmp(4,4,4) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 17. y^2 zx  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxyyz / (Sy(yy+tos)^2 * Sz(zz+tos) * Sx(xx+tos));
            val_tmp(2,2,2) = val_tmp(2,2,2) + tmp /4;
            val_tmp(4,2,2) = val_tmp(4,2,2) - tmp /4;
            val_tmp(2,2,4) = val_tmp(2,2,4) - tmp /4;
            val_tmp(4,2,4) = val_tmp(4,2,4) + tmp /4;
            val_tmp(2,3,2) = val_tmp(2,3,2) - tmp /2;
            val_tmp(4,3,2) = val_tmp(4,3,2) + tmp /2;
            val_tmp(2,3,4) = val_tmp(2,3,4) + tmp /2;
            val_tmp(4,3,4) = val_tmp(4,3,4) - tmp /2;
            val_tmp(2,4,2) = val_tmp(2,4,2) + tmp /4;
            val_tmp(4,4,2) = val_tmp(4,4,2) - tmp /4;
            val_tmp(2,4,4) = val_tmp(2,4,4) - tmp /4;
            val_tmp(4,4,4) = val_tmp(4,4,4) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 18. z^2 xy  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxyzz / (Sz(zz+tos)^2 * Sx(xx+tos) * Sy(yy+tos));
            val_tmp(2,2,2) = val_tmp(2,2,2) + tmp /4;
            val_tmp(4,2,2) = val_tmp(4,2,2) - tmp /4;
            val_tmp(2,4,2) = val_tmp(2,4,2) - tmp /4;
            val_tmp(4,4,2) = val_tmp(4,4,2) + tmp /4;
            val_tmp(2,2,3) = val_tmp(2,2,3) - tmp /2;
            val_tmp(4,2,3) = val_tmp(4,2,3) + tmp /2;
            val_tmp(2,4,3) = val_tmp(2,4,3) + tmp /2;
            val_tmp(4,4,3) = val_tmp(4,4,3) - tmp /2;
            val_tmp(2,2,4) = val_tmp(2,2,4) + tmp /4;
            val_tmp(4,2,4) = val_tmp(4,2,4) - tmp /4;
            val_tmp(2,4,4) = val_tmp(2,4,4) - tmp /4;
            val_tmp(4,4,4) = val_tmp(4,4,4) + tmp /4;
            
            
            %%%%%%%%%%%%%%%%%%%%% 19. x^2 y^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxyy / (Sx(xx+tos)*Sy(yy+tos))^2;
            val_tmp(2,2,3) = val_tmp(2,2,3) + tmp;
            val_tmp(3,2,3) = val_tmp(3,2,3) - tmp*2;
            val_tmp(4,2,3) = val_tmp(4,2,3) + tmp;
            val_tmp(2,3,3) = val_tmp(2,3,3) - tmp*2;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*4;
            val_tmp(4,3,3) = val_tmp(4,3,3) - tmp*2;
            val_tmp(2,4,3) = val_tmp(2,4,3) + tmp;
            val_tmp(3,4,3) = val_tmp(3,4,3) - tmp*2;
            val_tmp(4,4,3) = val_tmp(4,4,3) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%% 20. y^2 z^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cyyzz / (Sy(yy+tos)*Sz(zz+tos))^2;
            val_tmp(3,2,2) = val_tmp(3,2,2) + tmp;
            val_tmp(3,3,2) = val_tmp(3,3,2) - tmp*2;
            val_tmp(3,4,2) = val_tmp(3,4,2) + tmp;
            val_tmp(3,2,3) = val_tmp(3,2,3) - tmp*2;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*4;
            val_tmp(3,4,3) = val_tmp(3,4,3) - tmp*2;
            val_tmp(3,2,4) = val_tmp(3,2,4) + tmp;
            val_tmp(3,3,4) = val_tmp(3,3,4) - tmp*2;
            val_tmp(3,4,4) = val_tmp(3,4,4) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%% 21. z^2 x^2  %%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp = Cxxzz / (Sz(zz+tos)*Sx(xx+tos))^2;
            val_tmp(2,3,2) = val_tmp(2,3,2) + tmp;
            val_tmp(3,3,2) = val_tmp(3,3,2) - tmp*2;
            val_tmp(4,3,2) = val_tmp(4,3,2) + tmp;
            val_tmp(2,3,3) = val_tmp(2,3,3) - tmp*2;
            val_tmp(3,3,3) = val_tmp(3,3,3) + tmp*4;
            val_tmp(4,3,3) = val_tmp(4,3,3) - tmp*2;
            val_tmp(2,3,4) = val_tmp(2,3,4) + tmp;
            val_tmp(3,3,4) = val_tmp(3,3,4) - tmp*2;
            val_tmp(4,3,4) = val_tmp(4,3,4) + tmp;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%  mass term  %%%%%%%%%%%%%%%%%%%%%%%
            val_tmp(3,3,3) = val_tmp(3,3,3) - (hz*freq/vp)^2;
            

            for kk = 1:5
                for jj = 1:5
                    for ii = 1:5
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