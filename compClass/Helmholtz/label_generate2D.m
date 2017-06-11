function label = label_generate2D( Nz, Nx, TT, pts )

% zero out the label matrix
label = zeros(Nz,Nx);
counter = 0;

% Z is the fastest direction
for Ntmp = 1:length(TT)
    
    if mod(TT{Ntmp}.LL,2) == 1 && TT{Ntmp}.kids(1) ~= 0
        for zz = pts{Ntmp}.headZ : pts{Ntmp}.tailZ
            for xx = pts{Ntmp}.headX : pts{Ntmp}.tailX
                counter = counter + 1;
                label(zz,xx) = counter;
            end
        end
    else
        for xx = pts{Ntmp}.headX : pts{Ntmp}.tailX
            for zz = pts{Ntmp}.headZ : pts{Ntmp}.tailZ
                counter = counter + 1;
                label(zz,xx) = counter;
            end
        end
    end
    
end

return;
end