function pts = pts_generate3D( tos, Nx, Ny, Nz, TT, TL )

% sequence: X - Y - Z

Lvl = TT{1}.LL-1;
sample = struct('headX',0,'tailX',0,'headY',0,'tailY',0,'headZ',0,'tailZ',0);
pts = cell(1,length(TT));

% initialize
for Ntmp = 1:length(TT)
    pts{Ntmp} = struct(sample);
end
clear sample;

% initialize the root node
pts{end}.headX =  1;
pts{end}.tailX = Nx;
pts{end}.headY =  1;
pts{end}.tailY = Ny;
pts{end}.headZ =  1;
pts{end}.tailZ = Nz;

for level = 1:Lvl
    for ii = 1:length(TL{level})
        
        Ntmp = TL{level}(ii);
        kidL = TT{Ntmp}.kids(1);
        kidR = TT{Ntmp}.kids(2);
        
        if mod(level,3) == 1    % X direction
            
            pts{kidL}.headY = pts{Ntmp}.headY;
            pts{kidL}.tailY = pts{Ntmp}.tailY;
            pts{kidL}.headZ = pts{Ntmp}.headZ;
            pts{kidL}.tailZ = pts{Ntmp}.tailZ;
            pts{kidL}.headX = pts{Ntmp}.headX;
            
            pts{kidR}.headY = pts{Ntmp}.headY;
            pts{kidR}.tailY = pts{Ntmp}.tailY;
            pts{kidR}.headZ = pts{Ntmp}.headZ;
            pts{kidR}.tailZ = pts{Ntmp}.tailZ;
            pts{kidR}.tailX = pts{Ntmp}.tailX;
            
            pts{Ntmp}.headX = floor((pts{Ntmp}.headX + pts{Ntmp}.tailX - (tos-1))/2);
            pts{Ntmp}.tailX = pts{Ntmp}.headX + (tos-1);
            pts{kidL}.tailX = pts{Ntmp}.headX - 1;
            pts{kidR}.headX = pts{Ntmp}.tailX + 1;
            
        elseif mod(level,3) == 2    % Y direction
            
            pts{kidL}.headX = pts{Ntmp}.headX;
            pts{kidL}.tailX = pts{Ntmp}.tailX;
            pts{kidL}.headZ = pts{Ntmp}.headZ;
            pts{kidL}.tailZ = pts{Ntmp}.tailZ;
            pts{kidL}.headY = pts{Ntmp}.headY;
            
            pts{kidR}.headX = pts{Ntmp}.headX;
            pts{kidR}.tailX = pts{Ntmp}.tailX;
            pts{kidR}.headZ = pts{Ntmp}.headZ;
            pts{kidR}.tailZ = pts{Ntmp}.tailZ;
            pts{kidR}.tailY = pts{Ntmp}.tailY;
            
            pts{Ntmp}.headY = floor((pts{Ntmp}.headY + pts{Ntmp}.tailY - (tos-1))/2);
            pts{Ntmp}.tailY = pts{Ntmp}.headY + (tos-1);
            pts{kidL}.tailY = pts{Ntmp}.headY - 1;
            pts{kidR}.headY = pts{Ntmp}.tailY + 1;
            
        else    % Z direction
            
            pts{kidL}.headX = pts{Ntmp}.headX;
            pts{kidL}.tailX = pts{Ntmp}.tailX;
            pts{kidL}.headY = pts{Ntmp}.headY;
            pts{kidL}.tailY = pts{Ntmp}.tailY;
            pts{kidL}.headZ = pts{Ntmp}.headZ;
            
            pts{kidR}.headX = pts{Ntmp}.headX;
            pts{kidR}.tailX = pts{Ntmp}.tailX;
            pts{kidR}.headY = pts{Ntmp}.headY;
            pts{kidR}.tailY = pts{Ntmp}.tailY;
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