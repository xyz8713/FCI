function [ x, nflops ] = mat2hssParsolU( x, Q2, L, DL, U, B, V, hssTT, hssTL, rg )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  HSS parallel solution  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       A      *         x   =       b
% Q1' * A * Q2 * ( Q2' * x ) = Q1' * b

nflops = 0;

lenTT = length(hssTT);           % length of hssTT
Lvl   = length(hssTL);           % total level
SN    = lenTT - 1;               % Schur complement Node

%%%%%%%%%%  solution loop  %%%%%%%%%%
xt = cell(1,lenTT);   % for both RHS and solution

    
% INVERSE
for level = Lvl:-1:1
    for kk = 1:length(hssTL{level})

        Ntmp = hssTL{level}(kk);       % the current Node

        if hssTT{Ntmp}.kids(1) == 0    % base blocks where the data kicks in

            xt{Ntmp} = x( rg(1,Ntmp):rg(2,Ntmp) );   % data kicks in

        else

            kid1 = hssTT{Ntmp}.kids(1);
            kid2 = hssTT{Ntmp}.kids(2);
            rk1  = size( DL{kid1},1 );
            
            if kid2 == SN
                
                xt{Ntmp} = xt{kid1}( end-rk1+1:end );
                xt{Ntmp} = xt{Ntmp} - U{Ntmp} \ (L{Ntmp} \ (DL{Ntmp} * (U{kid1} * (B{kid1} * (V{kid2}' * xt{kid2})))));   % !!! IMPORTANT: square matrix
                nflops = nflops + flops('mv',U{kid1},'n') + flops('mv',B{kid1},'n') + flops('mv',V{kid2},'y');
                nflops = nflops + flops('sum', xt{Ntmp});
                break;
                
            end
            
            rk2  = size( DL{kid2},1 );
            xt{Ntmp} = [ xt{kid1}( end-rk1+1:end ); xt{kid2}( end-rk2+1:end ) ];

        end

    end
end


% FORWARD
for level = 1:Lvl
    for kk = 1:length(hssTL{level})

        Ntmp = hssTL{level}(kk);
        
        if Ntmp == lenTT
            
            rk1 = size( DL{kid1},1 );
            xt{kid1}( end-rk1+1:end ) = xt{Ntmp};
            
        elseif Ntmp == SN
            
            x( rg(1,Ntmp):rg(2,Ntmp) ) = xt{Ntmp};
            
        else
            
            xt{Ntmp} = Q2{Ntmp} * xt{Ntmp};
            nflops = nflops + flops('mv',Q2{Ntmp},'n');

            if hssTT{Ntmp}.kids(1) == 0
                x( rg(1,Ntmp):rg(2,Ntmp) ) = xt{Ntmp};
            else
                kid1 = hssTT{Ntmp}.kids(1);
                kid2 = hssTT{Ntmp}.kids(2);
                rk1 = size( DL{kid1},1 );
                rk2 = size( DL{kid2},1 );
                xt{kid1}( end-rk1+1:end ) = xt{Ntmp}( 1:rk1 );
                xt{kid2}( end-rk2+1:end ) = xt{Ntmp}( rk1+1:end );
            end
            
        end

    end
end


return;
end