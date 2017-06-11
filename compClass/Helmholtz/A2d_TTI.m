function  A = A2d_TTI( freq, tos, FS, Nz, Nx, N_delta, sigma0, hz, hx, label, vel_pz, rho, epsilon, delta, theta )
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
len = 25*Nz*Nx - 4*(16+26+9) - 30*((Nz-4)+(Nx-4));
row = zeros(1,len);
col = zeros(1,len);
val = zeros(1,len);
counter = 0;

% No use to follow the tree sequence
for xx = 1:Nx
    for zz = 1:Nz

        row_tmp = label(zz,xx)*ones(5,5);
        col_tmp = zeros(5,5);
        rho_tmp = ones(5,5);
        
        for jj = -2:2
            for ii = -2:2
                if zz+ii >= 1 && zz+ii <= Nz && xx+jj >= 1 && xx+jj <= Nx
                    col_tmp(ii+3,jj+3) = label(zz+ii,xx+jj);
                    rho_tmp(ii+3,jj+3) =   rho(zz+ii,xx+jj);
                end
            end
        end
               
        % value preprocessing
        vp  =  vel_pz(zz,xx);
        ee  = epsilon(zz,xx);
        dd  =   delta(zz,xx);
        tt  =   theta(zz,xx)*pi/180;
        D11 = cos(tt)^2;   D12 = sin(tt)^2;   D13 =  sin(2*tt);
        D21 = sin(tt)^2;   D22 = cos(tt)^2;   D23 = -sin(2*tt);
        
        Czz   =  -D11 - (1+2*ee)*D21;
        Cxx   = (-D12 - (1+2*ee)*D22) * rx^2;
        Czx   = (-D13 - (1+2*ee)*D23) * rx;
        Czzzz =  -2*(ee-dd)*(vp/hz/freq)^2*D11*D21;
        Cxxxx = (-2*(ee-dd)*(vp/hz/freq)^2*D12*D22) * rx^4;
        Czzzx = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D23+D21*D13)) * rx;
        Czxxx = (-2*(ee-dd)*(vp/hz/freq)^2*(D12*D23+D22*D13)) * rx^3;
        Czzxx = (-2*(ee-dd)*(vp/hz/freq)^2*(D11*D22+D21*D12+D13*D23)) * rx^2;
        
        % value generating
        val_tmp = zeros(5,5);
        
        %%%%%%%%%%%%%%%%%%%%%%  1. z^2  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czz * rho_tmp(3,3) / Sz(zz+tos);
        tmp1 = (1/(rho_tmp(2,3)*Sz(zz-1+tos)) + 1/(rho_tmp(3,3)*Sz(zz+tos)))/2;
        tmp2 = (1/(rho_tmp(4,3)*Sz(zz+1+tos)) + 1/(rho_tmp(3,3)*Sz(zz+tos)))/2;
        val_tmp(2,3) = val_tmp(2,3) + tmp * tmp1;
        val_tmp(3,3) = val_tmp(3,3) - tmp *(tmp1 + tmp2);
        val_tmp(4,3) = val_tmp(4,3) + tmp * tmp2;
        
        
        %%%%%%%%%%%%%%%%%%%%%%  2. x^2  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Cxx * rho_tmp(3,3) / Sx(xx+tos);
        tmp1 = (1/(rho_tmp(3,2)*Sx(xx-1+tos)) + 1/(rho_tmp(3,3)*Sx(xx+tos)))/2;
        tmp2 = (1/(rho_tmp(3,4)*Sx(xx+1+tos)) + 1/(rho_tmp(3,3)*Sx(xx+tos)))/2;
        val_tmp(3,2) = val_tmp(3,2) + tmp * tmp1;
        val_tmp(3,3) = val_tmp(3,3) - tmp *(tmp1 + tmp2);
        val_tmp(3,4) = val_tmp(3,4) + tmp * tmp2;
        

        %%%%%%%%%%%%%%%%%%%%%%  3. zx   %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czx * rho_tmp(3,3) / (4*Sz(zz+tos)*Sx(xx+tos));
        val_tmp(2,2) = val_tmp(2,2) + tmp * (rho_tmp(3,2)+rho_tmp(2,3))/(2*rho_tmp(3,2)*rho_tmp(2,3));
        val_tmp(4,2) = val_tmp(4,2) - tmp * (rho_tmp(3,2)+rho_tmp(4,3))/(2*rho_tmp(3,2)*rho_tmp(4,3));
        val_tmp(2,4) = val_tmp(2,4) - tmp * (rho_tmp(2,3)+rho_tmp(3,4))/(2*rho_tmp(2,3)*rho_tmp(3,4));
        val_tmp(4,4) = val_tmp(4,4) + tmp * (rho_tmp(4,3)+rho_tmp(3,4))/(2*rho_tmp(4,3)*rho_tmp(3,4));
        
        
        %%%%%%%%%%%%%%%%%%%%%%  4. z^4  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czzzz / Sz(zz+tos)^4;
        val_tmp(1,3) = val_tmp(1,3) + tmp;
        val_tmp(2,3) = val_tmp(2,3) - tmp*4;
        val_tmp(3,3) = val_tmp(3,3) + tmp*6;
        val_tmp(4,3) = val_tmp(4,3) - tmp*4;
        val_tmp(5,3) = val_tmp(5,3) + tmp;
        
        
        %%%%%%%%%%%%%%%%%%%%%%  5. x^4  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Cxxxx / Sx(xx+tos)^4;
        val_tmp(3,1) = val_tmp(3,1) + tmp;
        val_tmp(3,2) = val_tmp(3,2) - tmp*4;
        val_tmp(3,3) = val_tmp(3,3) + tmp*6;
        val_tmp(3,4) = val_tmp(3,4) - tmp*4;
        val_tmp(3,5) = val_tmp(3,5) + tmp;
        
        
        %%%%%%%%%%%%%%%%%%%%%%  6. z^3 x  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czzzx / (Sz(zz+tos)^3 * Sx(xx+tos));
        val_tmp(1,2) = val_tmp(1,2) + tmp / 4;
        val_tmp(2,2) = val_tmp(2,2) - tmp / 2;
        val_tmp(4,2) = val_tmp(4,2) + tmp / 2;
        val_tmp(5,2) = val_tmp(5,2) - tmp / 4;
        val_tmp(1,4) = val_tmp(1,4) - tmp / 4;
        val_tmp(2,4) = val_tmp(2,4) + tmp / 2;
        val_tmp(4,4) = val_tmp(4,4) - tmp / 2;
        val_tmp(5,4) = val_tmp(5,4) + tmp / 4;
        
        
        %%%%%%%%%%%%%%%%%%%%%%  7. z x^3  %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czxxx / (Sz(zz+tos) * Sx(xx+tos)^3);
        val_tmp(2,1) = val_tmp(2,1) + tmp / 4;
        val_tmp(2,2) = val_tmp(2,2) - tmp / 2;
        val_tmp(2,4) = val_tmp(2,4) + tmp / 2;
        val_tmp(2,5) = val_tmp(2,5) - tmp / 4;
        val_tmp(4,1) = val_tmp(4,1) - tmp / 4;
        val_tmp(4,2) = val_tmp(4,2) + tmp / 2;
        val_tmp(4,4) = val_tmp(4,4) - tmp / 2;
        val_tmp(4,5) = val_tmp(4,5) + tmp / 4;
        
        
        %%%%%%%%%%%%%%%%%%%%%%  8. z^2 x^2 %%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = Czzxx / (Sz(zz+tos)*Sx(xx+tos))^2;
        val_tmp(2,2) = val_tmp(2,2) + tmp;
        val_tmp(3,2) = val_tmp(3,2) - tmp*2;
        val_tmp(4,2) = val_tmp(4,2) + tmp;
        val_tmp(2,3) = val_tmp(2,3) - tmp*2;
        val_tmp(3,3) = val_tmp(3,3) + tmp*4;
        val_tmp(4,3) = val_tmp(4,3) - tmp*2;
        val_tmp(2,4) = val_tmp(2,4) + tmp;
        val_tmp(3,4) = val_tmp(3,4) - tmp*2;
        val_tmp(4,4) = val_tmp(4,4) + tmp;
        
        
        %%%%%%%%%%%%%%%%%%%%%  mass term  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        val_tmp(3,3) = val_tmp(3,3) - (hz*freq/vp)^2;
        

        for jj = 1:5
            for ii = 1:5
                if col_tmp(ii,jj) ~= 0
                    counter = counter+1;
                    row(counter) = row_tmp(ii,jj);
                    col(counter) = col_tmp(ii,jj);
                    val(counter) = val_tmp(ii,jj);
                end
            end
        end 
            
    end
end

A = sparse(row,col,val);

return;
end