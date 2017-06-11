function rg = rg_generate( m, hssTT )

% to calculate the range of each Node on the hssTT

rg = zeros( 2,length(hssTT) );

pt = 1;         % pointer of each basic block   
counter = 1;    % counter of basic blocks

for ii = 1:length(hssTT)
    
    if hssTT{ii}.kids(1) == 0 
        rg(1,ii) = pt;
        rg(2,ii) = pt + m(counter) - 1;
        pt = rg(2,ii) + 1; 
        counter = counter + 1;
    else
        rg(1,ii) = rg( 1, hssTT{ii}.kids(1) );
        rg(2,ii) = rg( 2, hssTT{ii}.kids(2) );
    end
    
end


return;
end