function label = label_generate3D( Nx, Ny, Nz, pts, TT, ST )

% label matrix
label = zeros(Nx, Ny, Nz);
counter = 0;

for Ntmp = 1:length(TT)
    
    if TT{Ntmp}.kids(1) == 0      % only for leaf nodes
        
        for zz = pts{Ntmp}.headZ : pts{Ntmp}.tailZ
            for yy = pts{Ntmp}.headY : pts{Ntmp}.tailY
                for xx = pts{Ntmp}.headX : pts{Ntmp}.tailX
                    counter = counter + 1;
                    label(xx,yy,zz) = counter;
                end
            end
        end
        
    else
        
        if mod(TT{Ntmp}.LL,3) == 1   % X direction
            
            for kk = 1:length(ST{Ntmp})
                node = ST{Ntmp}(kk);
                if mod(TT{node}.LL,3) == 0 && TT{node}.kids(1) ~= 0
                    for yy = pts{node}.headY : pts{node}.tailY
                        for zz = pts{node}.headZ : pts{node}.tailZ
                            for xx = pts{Ntmp}.headX : pts{Ntmp}.tailX
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                else
                    for zz = pts{node}.headZ : pts{node}.tailZ
                        for yy = pts{node}.headY : pts{node}.tailY
                            for xx = pts{Ntmp}.headX : pts{Ntmp}.tailX
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                end
            end
            
        elseif mod(TT{Ntmp}.LL,3) == 2   % Y direction
            
            for kk = 1:length(ST{Ntmp})
                node = ST{Ntmp}(kk);
                if mod(TT{node}.LL,3) == 0 && TT{node}.kids(1) ~= 0
                    for xx = pts{node}.headX : pts{node}.tailX
                        for zz = pts{node}.headZ : pts{node}.tailZ
                            for yy = pts{Ntmp}.headY : pts{Ntmp}.tailY
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                else
                    for zz = pts{node}.headZ : pts{node}.tailZ
                        for xx = pts{node}.headX : pts{node}.tailX
                            for yy = pts{Ntmp}.headY : pts{Ntmp}.tailY
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                end
            end
            
        else   % Z direction
            
            for kk = 1:length(ST{Ntmp})
                node = ST{Ntmp}(kk);
                if mod(TT{node}.LL,3) == 2 && TT{node}.kids(1) ~= 0
                    for xx = pts{node}.headX : pts{node}.tailX
                        for yy = pts{node}.headY : pts{node}.tailY
                            for zz = pts{Ntmp}.headZ : pts{Ntmp}.tailZ
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                else
                    for yy = pts{node}.headY : pts{node}.tailY
                        for xx = pts{node}.headX : pts{node}.tailX
                            for zz = pts{Ntmp}.headZ : pts{Ntmp}.tailZ
                                counter = counter + 1;
                                label(xx,yy,zz) = counter;
                            end
                        end
                    end
                end
            end
            
        end
        
    end
    
end


return;
end