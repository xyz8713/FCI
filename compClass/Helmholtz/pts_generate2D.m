function pts = pts_generate2D( tos, Nz, Nx, TT, TL )

% sequence: X - Z
% tos: Thickness of the Separator: 1 or 2

Lvl = TT{1}.LL-1;
sample = struct('headZ',0,'tailZ',0,'headX',0,'tailX',0);
pts = cell(1,length(TT));

% initialize
for Ntmp = 1:length(TT)
    pts{Ntmp} = struct(sample);
end
clear sample;

% initialize the root node
pts{end}.headZ =  1;
pts{end}.tailZ = Nz;
pts{end}.headX =  1;
pts{end}.tailX = Nx;

% domain decomposition loop by level
for level = 1:Lvl
    for ii = 1:length(TL{level})
        
        Ntmp = TL{level}(ii);
        kidL = TT{Ntmp}.kids(1);
        kidR = TT{Ntmp}.kids(2);
        
        if mod(level,2) == 1    % X direction
            
            pts{kidL}.headZ = pts{Ntmp}.headZ;
            pts{kidL}.tailZ = pts{Ntmp}.tailZ;
            pts{kidL}.headX = pts{Ntmp}.headX;
            
            pts{kidR}.headZ = pts{Ntmp}.headZ;
            pts{kidR}.tailZ = pts{Ntmp}.tailZ;
            pts{kidR}.tailX = pts{Ntmp}.tailX;
            
            pts{Ntmp}.headX = floor((pts{Ntmp}.headX + pts{Ntmp}.tailX - (tos-1))/2);
            pts{Ntmp}.tailX = pts{Ntmp}.headX + (tos-1);
            pts{kidL}.tailX = pts{Ntmp}.headX - 1;
            pts{kidR}.headX = pts{Ntmp}.tailX + 1;
            
        else    % Z direction
            
            pts{kidL}.headX = pts{Ntmp}.headX;
            pts{kidL}.tailX = pts{Ntmp}.tailX;
            pts{kidL}.headZ = pts{Ntmp}.headZ;
            
            pts{kidR}.headX = pts{Ntmp}.headX;
            pts{kidR}.tailX = pts{Ntmp}.tailX;
            pts{kidR}.tailZ = pts{Ntmp}.tailZ;
            
            pts{Ntmp}.headZ = floor((pts{Ntmp}.headZ + pts{Ntmp}.tailZ - (tos-1))/2);
            pts{Ntmp}.tailZ = pts{Ntmp}.headZ + (tos-1);
            pts{kidL}.tailZ = pts{Ntmp}.headZ - 1;
            pts{kidR}.headZ = pts{Ntmp}.tailZ + 1;
        
        end
        
    end
end


return;
end