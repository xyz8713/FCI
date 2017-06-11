function [ D, U, R, B, W, V, nflops ] = mat2hssParcmp( A, hssTT, hssTL, rg, tol )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  HSS parallel compression   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compr_ind = 0;    % 0: use matlab built-in function for QR with column pivoting
                  % 1: use compr_wang

nflops = 0;
lenTT = length(hssTT);     % length of the tree
Lvl   = hssTT{1}.LL;       % total level
D = cell(1,lenTT);
U = cell(1,lenTT);
R = cell(1,lenTT);
B = cell(1,lenTT);
W = cell(1,lenTT);
V = cell(1,lenTT);


%%%%%%%%%%%%%%%%%%%%%%  compression  %%%%%%%%%%%%%%%%%%%%%%%%

ls1 = zeros(1,lenTT);    % row blocks head point
ls2 = zeros(1,lenTT);    % row blocks tail point
for ii = 1:lenTT         % initialization
    ls1(ii) = rg(1,ii);
    ls2(ii) = rg(2,ii);
end

%%%%%%%%%%  1. row operation  %%%%%%%%%%
for level = Lvl:-1:2   % only to level 2: no need for level 1
    for kk = 1:length( hssTL{level} )
        
        Ntmp = hssTL{level}(kk);   % the current Node
        
        if hssTT{Ntmp}.kids(1) == 0   % only U factors
            
            %%%%%%%%%  DIAGONAL: D  %%%%%%%%%%
            D{Ntmp} = A( rg(1,Ntmp):rg(2,Ntmp), rg(1,Ntmp):rg(2,Ntmp) );
            
            % T1 are row blocks
            T1 = [ A( ls1(Ntmp):ls2(Ntmp), 1:rg(1,Ntmp)-1 ), A( ls1(Ntmp):ls2(Ntmp), rg(2,Ntmp)+1:end ) ];
            
            %%%%%%%%  COMPRESSION: U  %%%%%%%%%%
            [ U{Ntmp}, T1, nflops1 ] = compr_wang( T1, tol, compr_ind );   % overwrite T1
            nflops = nflops + nflops1; clear nflops1;
                
        else
            
            kid1 = hssTT{Ntmp}.kids(1);
            kid2 = hssTT{Ntmp}.kids(2);
            
            % move kid2 block upward to merge with kid1 block
            % all modifications are done on A
            ls2(Ntmp) = ls2(kid1) + ls2(kid2)-ls1(kid2) + 1;
            A( ls2(kid1)+1:ls2(Ntmp), 1:rg(1,Ntmp)-1 )   = A( ls1(kid2):ls2(kid2), 1:rg(1,Ntmp)-1 );
            A( ls2(kid1)+1:ls2(Ntmp), rg(2,Ntmp)+1:end ) = A( ls1(kid2):ls2(kid2), rg(2,Ntmp)+1:end );
            
            % merge horizontal two pieces
            T1 = [ A( ls1(Ntmp):ls2(Ntmp), 1:rg(1,Ntmp)-1 ), A( ls1(Ntmp):ls2(Ntmp), rg(2,Ntmp)+1:end ) ];
            
            %%%%%%%%  COMPRESSION: R  %%%%%%%%%%
            [ Rtmp, T1, nflops1 ] = compr_wang( T1, tol, compr_ind );   % overwrite T1
            nflops = nflops + nflops1; clear nflops1;
            
            R{kid1} = Rtmp( 1:ls2(kid1)-ls1(kid1)+1, : );
            R{kid2} = Rtmp( ls2(kid1)-ls1(kid1)+2:end, : );
            clear Rtmp;
            
        end
        
        % copy T1 back to A
        ls2(Ntmp) = ls1(Ntmp) + size(T1,1) - 1;
        A( ls1(Ntmp):ls2(Ntmp), 1:rg(1,Ntmp)-1 )   = T1( :, 1:rg(1,Ntmp)-1 );
        A( ls1(Ntmp):ls2(Ntmp), rg(2,Ntmp)+1:end ) = T1( :, rg(1,Ntmp):end );
        clear T1;
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ws1 = zeros(1,lenTT);
ws2 = zeros(1,lenTT);
for ii = 1:lenTT
    ws1(ii) = rg(1,ii);
    ws2(ii) = rg(2,ii);
end

%%%%%%%%%%  2. column operation  %%%%%%%%%%
for level = Lvl:-1:2
    for kk = 1:length( hssTL{level} )
        
        Ntmp = hssTL{level}(kk);
        
        %%%%  tmp is the collection of all neighbors  %%%%
        tmp = zeros(1,level-1);
        Ptmp = Ntmp;
        for ii = 1:length(tmp)
            if hssTT{Ptmp}.parent == Ptmp + 1  
                tmp(ii) = hssTT{ hssTT{Ptmp}.parent }.kids(1);
            else  
                tmp(ii) = hssTT{ hssTT{Ptmp}.parent }.kids(2);
            end
            Ptmp = hssTT{Ptmp}.parent;
        end
        tmp = sort(tmp,'ascend');  clear Ptmp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if hssTT{Ntmp}.kids(1) == 0   % only V factors
            
            counter = 1;
            for ii = 1:length(tmp)
                T2( counter:counter+ls2(tmp(ii))-ls1(tmp(ii)), 1:ws2(Ntmp)-ws1(Ntmp)+1 ) = A( ls1(tmp(ii)):ls2(tmp(ii)), ws1(Ntmp):ws2(Ntmp) );
                counter = counter + ls2(tmp(ii)) - ls1(tmp(ii)) + 1;
            end
            
            %%%%%%%%%  COMPRESSION: V  %%%%%%%%%%
            [ V{Ntmp}, T2, nflops1 ] = compr_wang( T2', tol, compr_ind );   % overwrite T2
            nflops = nflops + nflops1; clear nflops1;
            T2 = T2';
            
        else
            
            kid1 = hssTT{Ntmp}.kids(1);
            kid2 = hssTT{Ntmp}.kids(2);
            
            % move kid2 block leftward to merge with kid1
            ws2(Ntmp) = ws2(kid1) + ws2(kid2)-ws1(kid2) + 1;
            counter = 1;
            for ii = 1:length(tmp)
                A( ls1(tmp(ii)):ls2(tmp(ii)), ws2(kid1)+1:ws2(Ntmp) ) = A( ls1(tmp(ii)):ls2(tmp(ii)), ws1(kid2):ws2(kid2) );
                T2( counter:counter+ls2(tmp(ii))-ls1(tmp(ii)), 1:ws2(Ntmp)-ws1(Ntmp)+1 ) = A( ls1(tmp(ii)):ls2(tmp(ii)), ws1(Ntmp):ws2(Ntmp) );
                counter = counter + ls2(tmp(ii)) - ls1(tmp(ii)) + 1;
            end
            
            %%%%%%%%%%  COMPRESSION: W  %%%%%%%%%
            [ Wtmp, T2, nflops1 ] = compr_wang( T2', tol, compr_ind );   % overwrite T2
            nflops = nflops + nflops1; clear nflops1;
            T2 = T2';
            
            W{kid1} = Wtmp( 1:ws2(kid1)-ws1(kid1)+1, : );
            W{kid2} = Wtmp( ws2(kid1)-ws1(kid1)+2:end, : );
            clear Wtmp;
            
        end
            
        % copy T2 back to A
        ws2(Ntmp) = ws1(Ntmp) + size(T2,2) - 1;
        counter = 1;
        for ii = 1:length(tmp)
            if hssTT{tmp(ii)}.LL == level   % sibling Node
                %%%%%%%%%  PIVOTING NODE: B  %%%%%%%%%
                B{tmp(ii)} = T2( counter:counter+ls2(tmp(ii))-ls1(tmp(ii)), : );
            else
                A( ls1(tmp(ii)):ls2(tmp(ii)), ws1(Ntmp):ws2(Ntmp) ) = T2( counter:counter+ls2(tmp(ii))-ls1(tmp(ii)), : );
            end
            counter = counter + ls2(tmp(ii)) - ls1(tmp(ii)) + 1;
        end
        clear T2 tmp;
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
end