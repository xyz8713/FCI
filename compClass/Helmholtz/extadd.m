function [ UL, UR, LL, LR, nflops ] = extadd( UL, UR, LL, LR, Ptmp, U, Utmp )

% extend-add process for the current frontal matrix, which is partitioned
% into four individual parts:
% UL: upper left
% UR: upper right
% LL: lower left
% LR: lower right, which is the Schur complement

% size(UL) and length of the ind
sz  = Utmp{1}.Otail - Utmp{1}.Ohead + 1;
len = length(Utmp)-1;

% index of each block in U
ind = zeros(2,len);
counter = sz + 1;
for ii = 1:len
    ind(1,ii) = counter;
    ind(2,ii) = counter + Utmp{ii+1}.Otail - Utmp{ii+1}.Ohead;
    counter = ind(2,ii) + 1;
end

% mapping to the current frontal matrix
map = zeros(2,len);
for ii = 1:len
    counter = 1;
    for jj = 1:length(Ptmp)
        if Utmp{ii+1}.Ohead >= Ptmp{jj}.Ohead && Utmp{ii+1}.Otail <= Ptmp{jj}.Otail 
            map(1,ii) = Utmp{ii+1}.Ohead - Ptmp{jj}.Ohead + counter;
            map(2,ii) = map(1,ii) + ind(2,ii) - ind(1,ii);
            break;
        end
        counter = counter + Ptmp{jj}.Otail - Ptmp{jj}.Ohead + 1;
    end
end

% 1. extadd for UL
UL = UL + U(1:sz,1:sz);

% 2. extadd for UR
for jj = 1:len
    UR(1:sz,map(1,jj):map(2,jj)) = UR(1:sz,map(1,jj):map(2,jj)) + U(1:sz,ind(1,jj):ind(2,jj));
end

% 3. extadd for LL
for ii = 1:len
    LL(map(1,ii):map(2,ii),1:sz) = LL(map(1,ii):map(2,ii),1:sz) + U(ind(1,ii):ind(2,ii),1:sz);
end

% 4. extadd for LR
for jj = 1:len
    for ii = 1:len
        LR(map(1,ii):map(2,ii),map(1,jj):map(2,jj)) = LR(map(1,ii):map(2,ii),map(1,jj):map(2,jj)) + U(ind(1,ii):ind(2,ii),ind(1,jj):ind(2,jj));
        nflops = (ind(2,ii)-ind(1,ii)+1)*(ind(2,jj)-ind(1,jj));
    end
end


return;
end