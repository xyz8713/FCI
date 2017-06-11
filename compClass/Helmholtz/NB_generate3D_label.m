function nb = NB_generate3D_label( dim, tos, TT, pts, label )

% NB can be generated with label information
% sequence: X - Y- Z
nb = cell(1,length(TT));
sample = struct('Ohead',{0},'Otail',{0});     % overlapping head, overlapping tail, position in the matrix

for Ntmp = 1:length(TT)
    
    % search for the leftmost and rightmost kidL and kidR
    kidL = Ntmp;
    kidR = Ntmp;
    while TT{kidL}.kids(1) ~= 0
        kidL = TT{kidL}.kids(1);
        kidR = TT{kidR}.kids(2);
    end
    % the span of the current physical domain
    headXX = pts{kidL}.headX; headYY = pts{kidL}.headY; headZZ = pts{kidL}.headZ;
    tailXX = pts{kidR}.tailX; tailYY = pts{kidR}.tailY; tailZZ = pts{kidR}.tailZ;
    
    % calculating all neighbors
    ind = zeros(2,6);    % there are at most 6 neighbors
    counter = 0;
    Ptmp = TT{Ntmp}.parent;
    while Ptmp ~= 0
        if mod(TT{Ptmp}.LL,3) == 1    % X direction
            if headXX - pts{Ptmp}.tailX == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = -1;
            end
            if pts{Ptmp}.headX - tailXX == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = 1;
            end
        elseif mod(TT{Ptmp}.LL,3) == 2    % Y direction
            if headYY - pts{Ptmp}.tailY == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = -2;
            end
            if pts{Ptmp}.headY - tailYY == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = 2;
            end
        else    % Z direction
            if headZZ - pts{Ptmp}.tailZ == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = -3;
            end
            if pts{Ptmp}.headZ - tailZZ == 1
                counter = counter + 1;
                ind(1,counter) = Ptmp;
                ind(2,counter) = 3  ;
            end
        end
        Ptmp = TT{Ptmp}.parent;
    end
    
    % generating neighbors
    nb{Ntmp} = {};
    counter = 0;
    for ii = 1:6
        
        if ind(1,ii) == 0
            break;
        else
            Ptmp = ind(1,ii);
        end
        
        if abs(ind(2,ii)) == 1    % X direction
            
            counter = counter + 1;
            nb{Ntmp}{counter} = struct(sample);
            nb{Ntmp}{counter}.Ohead = (label( pts{Ptmp}.headX, headYY, headZZ ) - 1)*dim + 1;
            nb{Ntmp}{counter}.Otail = nb{Ntmp}{counter}.Ohead + (tailYY-headYY+1)*(tailZZ-headZZ+1)*tos*dim - 1;
            headY = headYY; tailY = tailYY; headZ = headZZ; tailZ = tailZZ;
            for jj = 1:ii
                if ind(2,jj) == -2
                    headY = headY - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( pts{Ptmp}.headX, headY, headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( pts{Ptmp}.tailX, headY+(tos-1), tailZ )*dim;
                end
                if ind(2,jj) == 2
                    tailY = tailY + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( pts{Ptmp}.headX, tailY-(tos-1), headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( pts{Ptmp}.tailX, tailY, tailZ )*dim;
                end
                if ind(2,jj) == -3
                    headZ = headZ - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( pts{Ptmp}.headX, headY, headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( pts{Ptmp}.tailX, tailY, headZ+(tos-1) )*dim;
                end
                if ind(2,jj) == 3
                    tailZ = tailZ + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( pts{Ptmp}.headX, headY, tailZ-(tos-1) )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( pts{Ptmp}.tailX, tailY, tailZ )*dim;
                end
            end
            
        elseif abs(ind(2,ii)) == 2    % Y direction
            
            counter = counter + 1;
            nb{Ntmp}{counter} = struct(sample);
            nb{Ntmp}{counter}.Ohead = (label( headXX, pts{Ptmp}.headY, headZZ )-1)*dim + 1;
            nb{Ntmp}{counter}.Otail = nb{Ntmp}{counter}.Ohead + (tailXX-headXX+1)*(tailZZ-headZZ+1)*tos*dim - 1;
            headX = headXX; tailX = tailXX; headZ = headZZ; tailZ = tailZZ;
            for jj = 1:ii
                if ind(2,jj) == -1
                    headX = headX - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, pts{Ptmp}.headY, headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( headX+(tos-1), pts{Ptmp}.tailY, tailZ )*dim;
                end
                if ind(2,jj) == 1
                    tailX = tailX + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( tailX-(tos-1), pts{Ptmp}.headY, headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, pts{Ptmp}.tailY, tailZ )*dim;
                end
                if ind(2,jj) == -3
                    headZ = headZ - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, pts{Ptmp}.headY, headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, pts{Ptmp}.tailY, headZ+(tos-1) )*dim;
                end
                if ind(2,jj) == 3
                    tailZ = tailZ + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, pts{Ptmp}.headY, tailZ-(tos-1) )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, pts{Ptmp}.tailY, tailZ )*dim;
                end
            end
            
        else    % Z direction
            
            counter = counter + 1;
            nb{Ntmp}{counter} = struct(sample);
            nb{Ntmp}{counter}.Ohead = (label( headXX, headYY, pts{Ptmp}.headZ )-1)*dim + 1;
            nb{Ntmp}{counter}.Otail = nb{Ntmp}{counter}.Ohead + (tailXX-headXX+1)*(tailYY-headYY+1)*tos*dim - 1;
            headX = headXX; tailX = tailXX; headY = headYY; tailY = tailYY;
            for jj = 1:ii
                if ind(2,jj) == -1
                    headX = headX - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, headY, pts{Ptmp}.headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( headX+(tos-1), tailY, pts{Ptmp}.tailZ )*dim;
                end
                if ind(2,jj) == 1
                    tailX = tailX + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( tailX-(tos-1), headY, pts{Ptmp}.headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, tailY, pts{Ptmp}.tailZ )*dim;
                end
                if ind(2,jj) == -2
                    headY = headY - tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, headY, pts{Ptmp}.headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, headY+(tos-1), pts{Ptmp}.tailZ )*dim;
                end
                if ind(2,jj) == 2
                    tailY = tailY + tos;
                    counter = counter + 1;
                    nb{Ntmp}{counter} = struct(sample);
                    nb{Ntmp}{counter}.Ohead = (label( headX, tailY-(tos-1), pts{Ptmp}.headZ )-1)*dim + 1;
                    nb{Ntmp}{counter}.Otail = label( tailX, tailY, pts{Ptmp}.tailZ )*dim;
                end
            end
            
        end 
        
    end
        
end


return;
end