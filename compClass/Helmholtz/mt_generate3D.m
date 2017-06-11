function mt = mt_generate3D( TT, Tlabel, NB )

% mt is the fraction of the leading block in the process of multifrontal
Lvl = TT{1}.LL-1;
mt = cell(1,length(TT));

for Ntmp = 1:length(TT)
    
    counter = 0;
    level = TT{Ntmp}.LL + 1;
    while level <= Lvl
        if mod(level-TT{Ntmp}.LL,3) ~= 0
            counter = counter + 1;
        end
        level = level + 1;
    end
    Nfrac = 2^counter;    % number of fractions
    
    range = Tlabel(Ntmp+1) - Tlabel(Ntmp);
    frac = floor(range/Nfrac);
    
    if Ntmp ~= length(TT)
        
        mt{Ntmp} = frac*ones(1,Nfrac+1);
        mt{Ntmp}(Nfrac) = range - frac*(Nfrac-1);
        counter = 0;  % Schur complement is the last block
        for ii = 1:length(NB{Ntmp})
            counter = counter + ( NB{Ntmp}{ii}.Otail - NB{Ntmp}{ii}.Ohead + 1 );    % only the overlapping parts account for
        end
        mt{Ntmp}(Nfrac+1) = counter;
        
    else
        
        mt{Ntmp} = frac*ones(1,Nfrac);   % NO Schur complement
        mt{Ntmp}(Nfrac) = range - frac*(Nfrac-1);
        
    end;
    
end


return;
end