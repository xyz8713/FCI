function index = index3D( xx, yy, zz, Nx, Ny, Nz, TT, ST, pts, Tlabel )


if (xx<2) || (xx>Nx-1) || (yy<2) || (yy>Ny-1) || (zz<2) || (zz>Nz-1)
    
    index = 0; return;
    
else

    Ntmp = length(TT);   % starting from the root Node

    % 3D binary tree search
    while TT{Ntmp}.kids(1) ~= 0

        if mod(TT{Ntmp}.LL,3) == 1    % X direction

            if xx < pts{Ntmp}.headX
                Ntmp = TT{Ntmp}.kids(1);
            elseif xx > pts{Ntmp}.headX
                Ntmp = TT{Ntmp}.kids(2);
            else 
                
                % 2D tree preprocessing
                root = length(ST{Ntmp});      % 2D subtree imbedded
                lvl  = TT{root}.LL-1;
                Tlabel2D = zeros(1,root+1);   % 2D sub Tlabel
                Tlabel2D(1) = Tlabel(Ntmp);
                for kk = 1:root
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ST{Ntmp}(kk)}.tailY-pts{ST{Ntmp}(kk)}.headY+1)*(pts{ST{Ntmp}(kk)}.tailZ-pts{ST{Ntmp}(kk)}.headZ+1);
                end

                % 2D binary tree search
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    if mod(TT{kk}.LL-lvl,2) == 1   % Y direction
                        if yy < pts{ST{Ntmp}(kk)}.headY
                            kk = TT{kk}.kids(1);
                        elseif yy > pts{ST{Ntmp}(kk)}.headY
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + zz - pts{ST{Ntmp}(kk)}.headZ;
                            return;
                        end
                    else    % Z direction
                        if zz < pts{ST{Ntmp}(kk)}.headZ
                            kk = TT{kk}.kids(1);
                        elseif zz > pts{ST{Ntmp}(kk)}.headZ
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + yy - pts{ST{Ntmp}(kk)}.headY;
                            return;
                        end
                    end
                end
                index = Tlabel2D(kk) + (zz-pts{ST{Ntmp}(kk)}.headZ)*(pts{ST{Ntmp}(kk)}.tailY-pts{ST{Ntmp}(kk)}.headY+1) + (yy-pts{ST{Ntmp}(kk)}.headY);
                return;
                
            end

        elseif mod(TT{Ntmp}.LL,3) == 2    % Y direction
            
            if yy < pts{Ntmp}.headY
                Ntmp = TT{Ntmp}.kids(1);
            elseif yy > pts{Ntmp}.headY
                Ntmp = TT{Ntmp}.kids(2);
            else
                
                % 2D tree preprocessing
                root = length(ST{Ntmp});      % 2D subtree imbedded
                lvl  = TT{root}.LL-1;
                Tlabel2D = zeros(1,root+1);   % 2D sub Tlabel
                Tlabel2D(1) = Tlabel(Ntmp);
                for kk = 1:root
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ST{Ntmp}(kk)}.tailX-pts{ST{Ntmp}(kk)}.headX+1)*(pts{ST{Ntmp}(kk)}.tailZ-pts{ST{Ntmp}(kk)}.headZ+1);
                end

                % 2D binary tree search
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    if mod(TT{kk}.LL-lvl,2) == 1   % Z direction
                        if zz < pts{ST{Ntmp}(kk)}.headZ
                            kk = TT{kk}.kids(1);
                        elseif zz > pts{ST{Ntmp}(kk)}.headZ
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + xx - pts{ST{Ntmp}(kk)}.headX;
                            return;
                        end
                    else    % X direction
                        if xx < pts{ST{Ntmp}(kk)}.headX
                            kk = TT{kk}.kids(1);
                        elseif xx > pts{ST{Ntmp}(kk)}.headX
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + zz - pts{ST{Ntmp}(kk)}.headZ;
                            return;
                        end
                    end
                end
                index = Tlabel2D(kk) + (zz-pts{ST{Ntmp}(kk)}.headZ)*(pts{ST{Ntmp}(kk)}.tailX-pts{ST{Ntmp}(kk)}.headX+1) + (xx-pts{ST{Ntmp}(kk)}.headX);
                return;
                
            end

        else    % Z direction

            if zz < pts{Ntmp}.headZ
                Ntmp = TT{Ntmp}.kids(1);
            elseif zz > pts{Ntmp}.headZ
                Ntmp = TT{Ntmp}.kids(2);
            else
                
                % 2D tree preprocessing
                root = length(ST{Ntmp});      % 2D subtree imbedded
                lvl  = TT{root}.LL-1;
                Tlabel2D = zeros(1,root+1);   % 2D sub Tlabel
                Tlabel2D(1) = Tlabel(Ntmp);
                for kk = 1:root
                    Tlabel2D(kk+1) = Tlabel2D(kk) + (pts{ST{Ntmp}(kk)}.tailX-pts{ST{Ntmp}(kk)}.headX+1)*(pts{ST{Ntmp}(kk)}.tailY-pts{ST{Ntmp}(kk)}.headY+1);
                end

                % 2D binary tree search
                kk = root;
                while TT{kk}.kids(1) ~= 0
                    if mod(TT{kk}.LL-lvl,2) == 1   % X direction
                        if xx < pts{ST{Ntmp}(kk)}.headX
                            kk = TT{kk}.kids(1);
                        elseif xx > pts{ST{Ntmp}(kk)}.headX
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + yy - pts{ST{Ntmp}(kk)}.headY;
                            return;
                        end
                    else    % Y direction
                        if yy < pts{ST{Ntmp}(kk)}.headY
                            kk = TT{kk}.kids(1);
                        elseif yy > pts{ST{Ntmp}(kk)}.headY
                            kk = TT{kk}.kids(2);
                        else
                            index = Tlabel2D(kk) + xx - pts{ST{Ntmp}(kk)}.headX;
                            return;
                        end
                    end
                end
                index = Tlabel2D(kk) + (yy-pts{ST{Ntmp}(kk)}.headY)*(pts{ST{Ntmp}(kk)}.tailX-pts{ST{Ntmp}(kk)}.headX+1) + (xx-pts{ST{Ntmp}(kk)}.headX);
                return;
                
            end

        end

    end

    index = Tlabel(Ntmp) + (zz-pts{Ntmp}.headZ)*(pts{Ntmp}.tailX-pts{Ntmp}.headX+1)*(pts{Ntmp}.tailY-pts{Ntmp}.headY+1) ...
                         + (yy-pts{Ntmp}.headY)*(pts{Ntmp}.tailX-pts{Ntmp}.headX+1) ...
                         + (xx-pts{Ntmp}.headX);
    return;

end


end