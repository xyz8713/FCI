function nb = NB_generate3D( dim, tos, TT, pts, Tlabel, ST )

% NB can be generated without knowing label information
% sequence: X - Y - Z
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
    
    nb{Ntmp} = {};    % a doubly linked list
    counter = 0;
    Ptmp = TT{Ntmp}.parent;
    
    while Ptmp ~= 0
        
        if mod(TT{Ptmp}.LL,3) == 1    % X direction
            if (headXX - pts{Ptmp}.tailX == 1) || (pts{Ptmp}.headX - tailXX == 1)
                
                % 1. 2D tree preprocessing
                root = length(ST{Ptmp});      % 2D subtree is imbedded in the global tree
                lvl  = TT{root}.LL-1;         % level waiting to be subtracted
                Tlabel2D = zeros(1,root+1);   % 2D subTlabel
                Tlabel2D(1) = Tlabel(Ptmp);
                for kk = 1:root
                    ind = ST{Ptmp}(kk);
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ind}.tailY-pts{ind}.headY+1)*(pts{ind}.tailZ-pts{ind}.headZ+1)*tos*dim;
                end
                
                % 2. search for the node which spans the 2D projected domain
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    ind = ST{Ptmp}(kk);
                    if mod(TT{kk}.LL-lvl,2) == 1    % Y direction
                        if pts{ind}.headY > tailYY
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailY < headYY
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    else    % Z direction
                        if pts{ind}.headZ > tailZZ
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailZ < headZZ
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    end
                end
                
                % 3. search for the leftmost kid
                kidL = kk;
                while TT{kidL}.kids(1) ~= 0
                    kidL = TT{kidL}.kids(1);
                end
                counter = counter + 1;
                nb{Ntmp}{counter} = struct(sample);
                nb{Ntmp}{counter}.Ohead = Tlabel2D(kidL);
                nb{Ntmp}{counter}.Otail = Tlabel2D(kk+1) - 1;
                
                % 4. 2D search for neighbors
                P2D = TT{kk}.parent;
                headY = headYY; tailY = tailYY; headZ = headZZ; tailZ = tailZZ;
                
                while P2D ~= TT{root}.parent
                    
                    ind = ST{Ptmp}(P2D);
                    
                    if mod(TT{P2D}.LL-lvl,2) == 1   % Y direction
                        
                        if headY - pts{ind}.tailY == 1
                            headY = headY - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headZ  -pts{ind}.headZ)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailZ+1-pts{ind}.headZ)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headY - tailY == 1
                            tailY = tailY + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headZ  -pts{ind}.headZ)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailZ+1-pts{ind}.headZ)*tos^2*dim - 1;
                        end
                        
                    else    % Z direction
                        
                        if headZ - pts{ind}.tailZ == 1
                            headZ = headZ - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headY  -pts{ind}.headY)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailY+1-pts{ind}.headY)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headZ - tailZ == 1
                            tailZ = tailZ + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headY  -pts{ind}.headY)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailY+1-pts{ind}.headY)*tos^2*dim - 1;
                        end
                        
                    end
                    
                    P2D = TT{P2D}.parent;
                    
                end
                clear Tlabel2D;
            end
                
            
        elseif mod(TT{Ptmp}.LL,3) == 2    % Y direction
            if (headYY - pts{Ptmp}.tailY == 1) || (pts{Ptmp}.headY - tailYY == 1)
                
                % 1. 2D tree preprocessing
                root = length(ST{Ptmp});      % 2D subtree is imbedded in the global tree
                lvl  = TT{root}.LL-1;         % level waiting to be subtracted
                Tlabel2D = zeros(1,root+1);   % 2D subTlabel
                Tlabel2D(1) = Tlabel(Ptmp);
                for kk = 1:root
                    ind = ST{Ptmp}(kk);
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ind}.tailX-pts{ind}.headX+1)*(pts{ind}.tailZ-pts{ind}.headZ+1)*tos*dim;
                end
                
                % 2. search for the node which spans the 2D projected domain
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    ind = ST{Ptmp}(kk);
                    if mod(TT{kk}.LL-lvl,2) == 1    % Z direction
                        if pts{ind}.headZ > tailZZ
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailZ < headZZ
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    else    % X direction
                        if pts{ind}.headX > tailXX
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailX < headXX
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    end
                end
                
                % 3. search for the leftmost kid
                kidL = kk;
                while TT{kidL}.kids(1) ~= 0
                    kidL = TT{kidL}.kids(1);
                end
                counter = counter + 1;
                nb{Ntmp}{counter} = struct(sample);
                nb{Ntmp}{counter}.Ohead = Tlabel2D(kidL);
                nb{Ntmp}{counter}.Otail = Tlabel2D(kk+1) - 1;
                
                % 4. 2D search for neighbors
                P2D = TT{kk}.parent;
                headX = headXX; tailX = tailXX; headZ = headZZ; tailZ = tailZZ;
                
                while P2D ~= TT{root}.parent
                    
                    ind = ST{Ptmp}(P2D);
                    
                    if mod(TT{P2D}.LL-lvl,2) == 1   % Z direction
                        
                        if headZ - pts{ind}.tailZ == 1
                            headZ = headZ - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headX  -pts{ind}.headX)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailX+1-pts{ind}.headX)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headZ - tailZ == 1
                            tailZ = tailZ + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headX  -pts{ind}.headX)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailX+1-pts{ind}.headX)*tos^2*dim - 1;
                        end
                        
                    else    % X direction
                        
                        if headX - pts{ind}.tailX == 1
                            headX = headX - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headZ  -pts{ind}.headZ)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailZ+1-pts{ind}.headZ)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headX - tailX == 1
                            tailX = tailX + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headZ  -pts{ind}.headZ)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailZ+1-pts{ind}.headZ)*tos^2*dim - 1;
                        end
                        
                    end
                    
                    P2D = TT{P2D}.parent;
                    
                end
                clear Tlabel2D;
            end
            
        else    % Z direction
            if (headZZ - pts{Ptmp}.tailZ == 1) || (pts{Ptmp}.headZ - tailZZ == 1)
                
                 % 1. 2D tree preprocessing
                root = length(ST{Ptmp});      % 2D subtree is imbedded in the global tree
                lvl  = TT{root}.LL-1;         % level waiting to be subtracted
                Tlabel2D = zeros(1,root+1);   % 2D subTlabel
                Tlabel2D(1) = Tlabel(Ptmp);
                for kk = 1:root
                    ind = ST{Ptmp}(kk);
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ind}.tailX-pts{ind}.headX+1)*(pts{ind}.tailY-pts{ind}.headY+1)*tos*dim;
                end
                
                % 2. search for the node which spans the 2D projected domain
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    ind = ST{Ptmp}(kk);
                    if mod(TT{kk}.LL-lvl,2) == 1    % X direction
                        if pts{ind}.headX > tailXX
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailX < headXX
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    else    % Y direction
                        if pts{ind}.headY > tailYY
                            kk = TT{kk}.kids(1);
                        elseif pts{ind}.tailY < headYY
                            kk = TT{kk}.kids(2);
                        else
                            break;
                        end
                    end
                end
                
                % 3. search for the leftmost kid
                kidL = kk;
                while TT{kidL}.kids(1) ~= 0
                    kidL = TT{kidL}.kids(1);
                end
                counter = counter + 1;
                nb{Ntmp}{counter} = struct(sample);
                nb{Ntmp}{counter}.Ohead = Tlabel2D(kidL);
                nb{Ntmp}{counter}.Otail = Tlabel2D(kk+1) - 1;
                
                % 4. 2D search for neighbors
                P2D = TT{kk}.parent;
                headX = headXX; tailX = tailXX; headY = headYY; tailY = tailYY;
                
                while P2D ~= TT{root}.parent
                    
                    ind = ST{Ptmp}(P2D);
                    
                    if mod(TT{P2D}.LL-lvl,2) == 1   % X direction
                        
                        if headX - pts{ind}.tailX == 1
                            headX = headX - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headY  -pts{ind}.headY)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailY+1-pts{ind}.headY)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headX - tailX == 1
                            tailX = tailX + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headY  -pts{ind}.headY)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailY+1-pts{ind}.headY)*tos^2*dim - 1;
                        end
                        
                    else    % Y direction
                        
                        if headY - pts{ind}.tailY == 1
                            headY = headY - tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headX  -pts{ind}.headX)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailX+1-pts{ind}.headX)*tos^2*dim - 1;
                        end
                        
                        if pts{ind}.headY - tailY == 1
                            tailY = tailY + tos;
                            counter = counter + 1;
                            nb{Ntmp}{counter} = struct(sample);
                            nb{Ntmp}{counter}.Ohead = Tlabel2D(P2D) + (headX  -pts{ind}.headX)*tos^2*dim;
                            nb{Ntmp}{counter}.Otail = Tlabel2D(P2D) + (tailX+1-pts{ind}.headX)*tos^2*dim - 1;
                        end
                        
                    end
                    
                    P2D = TT{P2D}.parent;
                    
                end
                clear Tlabel2D;
            end
            
        end
        
        Ptmp = TT{Ptmp}.parent;
        
    end
    
end


return;
end