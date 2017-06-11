function ST = ST_generate3D( TT, pts )

% SubTree generation for 3D case: each 2D separator contains a
% substructure
Lvl = TT{1}.LL;
ST = cell(1,length(TT));

for Ntmp = 1:length(TT)
    if TT{Ntmp}.kids(1) ~= 0
        
        % computing the total level of the SubTree first
        counter = 0;
        for level = TT{Ntmp}.LL+1:Lvl
            if mod(level-TT{Ntmp}.LL,3) == 0
                if level == Lvl
                    counter = counter + 1;
                    break;
                else
                    continue;
                end
            else
                counter = counter + 1;
            end
        end
        ST{Ntmp} = zeros(1,2^counter-1);   % allocate space
        
        % search for the leftmost kid
        kidL = TT{Ntmp}.kids(1);
        while TT{kidL}.kids(1) ~= 0
            kidL = TT{kidL}.kids(1);
        end
        
        % build up the SubTree
        counter = 0;
        if mod(TT{Ntmp}.LL,3) == 1    % X direction separators
            for kk = kidL:TT{Ntmp}.kids(1)
                if pts{Ntmp}.headX - pts{kk}.tailX == 1
                    counter = counter + 1;
                    ST{Ntmp}(counter) = kk;
                end
            end
        elseif mod(TT{Ntmp}.LL,3) == 2    % Y direction separators
            for kk = kidL:TT{Ntmp}.kids(1)
                if pts{Ntmp}.headY - pts{kk}.tailY == 1
                    counter = counter + 1;
                    ST{Ntmp}(counter) = kk;
                end
            end
        else    % Z direction separators
            for kk = kidL:TT{Ntmp}.kids(1)
                if pts{Ntmp}.headZ - pts{kk}.tailZ == 1
                    counter = counter + 1;
                    ST{Ntmp}(counter) = kk;
                end
            end
        end
        
    end
end


return;
end