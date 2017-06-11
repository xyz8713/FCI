function Tlabel = Tlabel_generate3D( dim, pts )

% records the starting and ending points of each Node in the full matrix A

Tlabel = zeros(1,length(pts)+1);
Tlabel(1) = 1;

for Ntmp = 1:length(pts)   % this sequence is important
    Tlabel(Ntmp+1) = Tlabel(Ntmp) + (pts{Ntmp}.tailX-pts{Ntmp}.headX+1) * (pts{Ntmp}.tailY-pts{Ntmp}.headY+1) * (pts{Ntmp}.tailZ-pts{Ntmp}.headZ+1) * dim;
end

return;
end