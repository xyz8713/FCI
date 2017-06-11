function nb = NB_generate2D( dim, tos, TT, pts, Tlabel )

% tos: Thickness of Separator

nb = cell(1,length(TT));
sample = struct('Ohead',{0},'Otail',{0});     % overlapping head, overlapping tail, position in the matrix

for Ntmp = 1:length(TT)
    
    % search for leftmost and rightmost kids: kidL and kidR span the whole physical domain
    kidL = Ntmp;    % leftmost  kid 
    kidR = Ntmp;    % rightmost kid
    while TT{kidL}.kids(1) ~= 0
        kidL = TT{kidL}.kids(1);
        kidR = TT{kidR}.kids(2);
    end
    headZZ = pts{kidL}.headZ;  
    headXX = pts{kidL}.headX;
    tailZZ = pts{kidR}.tailZ;  
    tailXX = pts{kidR}.tailX;
    
    nb{Ntmp} = {};   % a doubly linked list
    ii = 0;   % counter for neighbors
    Ptmp = TT{Ntmp}.parent;  % initialize the current parent neighbor
    
    while Ptmp ~= 0
        
        if mod(TT{Ptmp}.LL,2) == 1    % X direction separators
            
            if headXX - pts{Ptmp}.tailX == 1
                ii = ii + 1;
                headXX = headXX - tos;
                nb{Ntmp}{ii} = struct(sample);
                nb{Ntmp}{ii}.Ohead = Tlabel(Ptmp) + (headZZ  -pts{Ptmp}.headZ)*tos*dim;
                nb{Ntmp}{ii}.Otail = Tlabel(Ptmp) + (tailZZ+1-pts{Ptmp}.headZ)*tos*dim - 1;
            end

            if pts{Ptmp}.headX - tailXX == 1
                ii = ii + 1;
                tailXX = tailXX + tos;
                nb{Ntmp}{ii} = struct(sample);
                nb{Ntmp}{ii}.Ohead = Tlabel(Ptmp) + (headZZ  -pts{Ptmp}.headZ)*tos*dim;
                nb{Ntmp}{ii}.Otail = Tlabel(Ptmp) + (tailZZ+1-pts{Ptmp}.headZ)*tos*dim - 1;
            end
            
        else    % Z direction         
        
            if headZZ - pts{Ptmp}.tailZ == 1
                ii = ii + 1;
                headZZ = headZZ - tos; 
                nb{Ntmp}{ii} = struct(sample);
                nb{Ntmp}{ii}.Ohead = Tlabel(Ptmp) + (headXX  -pts{Ptmp}.headX)*tos*dim;
                nb{Ntmp}{ii}.Otail = Tlabel(Ptmp) + (tailXX+1-pts{Ptmp}.headX)*tos*dim - 1;
            end

            if pts{Ptmp}.headZ - tailZZ == 1
                ii = ii + 1;
                tailZZ = tailZZ + tos;
                nb{Ntmp}{ii} = struct(sample);
                nb{Ntmp}{ii}.Ohead = Tlabel(Ptmp) + (headXX  -pts{Ptmp}.headX)*tos*dim;
                nb{Ntmp}{ii}.Otail = Tlabel(Ptmp) + (tailXX+1-pts{Ptmp}.headX)*tos*dim - 1;
            end
            
        end
            
        Ptmp = TT{Ptmp}.parent;
        
    end
    
end


return;
end