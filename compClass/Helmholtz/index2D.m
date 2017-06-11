function index = index2D( zz, xx, Nz, Nx, TT, pts, Tlabel )


if (zz<2) || (zz>Nz-1) || (xx<2) || (xx>Nx-1)    % outside the computational domain
    
    index = 0; return;
    
else

    Ntmp = length(TT);   % starting from the root Node

    % 2D binary tree search
    while TT{Ntmp}.kids(1) ~= 0

        if mod(TT{Ntmp}.LL,2) == 1   % X direction
            if xx < pts{Ntmp}.headX
                Ntmp = TT{Ntmp}.kids(1);
            elseif xx > pts{Ntmp}.headX
                Ntmp = TT{Ntmp}.kids(2);
            else
                index = Tlabel(Ntmp) + zz - pts{Ntmp}.headZ;
                return;
            end
        else   % Z direction
            if zz < pts{Ntmp}.headZ
                Ntmp = TT{Ntmp}.kids(1);
            elseif zz > pts{Ntmp}.headZ
                Ntmp = TT{Ntmp}.kids(2);
            else
                index = Tlabel(Ntmp) + xx - pts{Ntmp}.headX;
                return;
            end
        end

    end
    
    index = Tlabel(Ntmp) + (xx-pts{Ntmp}.headX)*(pts{Ntmp}.tailZ-pts{Ntmp}.headZ+1) + (zz-pts{Ntmp}.headZ);
    return;
    
end


end