function TL = TL_generate(TT)

TL = cell( 1, TT{1}.LL );
counter = zeros(1,length(TL));

for Ntmp = 1:length(TT)
    
    tmp = TT{Ntmp}.LL;
    counter(tmp) = counter(tmp) + 1;
    TL{tmp}( counter(tmp) ) = Ntmp;
    
end

return;
end