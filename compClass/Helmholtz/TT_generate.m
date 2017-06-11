function TT = TT_generate( N )

% to generate the postordering binary tree based on the given number of
% total nodes N

% 1. upper level bound
L1 = floor( log(N+1)/log(2) );   % level lower bound
L2 = ceil( log(N+1)/log(2) );    % level upper bound

% 2. initialize TT
sample = struct('LL',{0},'parent',{0},'kids',{[0,0]});
TT = cell(1,N);  
for kk = 1:N
    TT{kk} = struct(sample);
end
clear sample;

% 3. build up TT: without considering any level info
nop = zeros(L2,3);
for kk = 1:N   % loop over all nodes
    ind = 0;
    for level = 1:L2
        if nop(level,1) == 2
            ind = level;        % which level has two no-parent nodes
            break;
        end
    end
    
    if ind == 0              % go the the lowest level
        level = L2;
        nop(level,1) = nop(level,1) + 1;
        if nop(level,2) == 0
            nop(level,2) = kk;
        else
            nop(level,3) = kk;
        end
    else
        level = ind - 1;     % go to the upper level
        nop(level,1) = nop(level,1) + 1;
        if nop(level,2) == 0
            nop(level,2) = kk;
        else
            nop(level,3) = kk;
        end
        TT{kk}.kids = nop(ind,2:3);
        TT{nop(ind,2)}.parent = kk;
        TT{nop(ind,3)}.parent = kk;
        nop(ind,:) = zeros(1,3);      % clear step
    end
    
end
clear nop;

% 4. modify TT connections
if L1 ~= L2
    
    % initialize tmp: calculate nodes with ( parent == 0 )
    tmp = [];
    counter = 0;
    for kk = 1:N
        if TT{kk}.parent == 0
            counter = counter + 1;
            tmp(counter) = kk;
        end
    end
    
    while length(tmp) > 1
        
        Node = tmp(end);  % out the stack
        if TT{Node}.kids(1) ~= 0   % in the stack
            tmp(end)   = TT{Node}.kids(1);
            tmp(end+1) = TT{Node}.kids(2);
        else
            tmp = tmp(1:end-1);
        end
        TT{Node}.kids(1) = tmp(1);
        TT{Node}.kids(2) = tmp(end);
        TT{ TT{Node}.kids(1) }.parent = Node;
        TT{ TT{Node}.kids(2) }.parent = Node;
        tmp = tmp(2:end);
        
    end
    
    clear tmp Node counter;
    
end

% 5. set up level info for the TT
TT{N}.LL = 1;
for kk = N-1:-1:1
    TT{kk}.LL = TT{ TT{kk}.parent }.LL + 1;
end


return;
end