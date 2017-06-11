function [ x, nflops ] = mat2hssParsolL( x, Q1, L, DL, U, B, W, V, hssTT, hssTL, rg )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  HSS parallel solution: L part  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       A      *         x   =       b
% Q1' * A * Q2 * ( Q2' * x ) = Q1' * b

nflops = 0;

lenTT = length(hssTT);           % length of hssTT
Lvl   = length(hssTL);           % total level
SN    = lenTT - 1;               % Schur complement Node

%%%%%%%%%%  solution loop  %%%%%%%%%%
xt = cell(1,lenTT);   % for both RHS and solution
Vx = cell(1,lenTT);   % MOST IMPORTANT OBJECT!!!!!!!

    
% INVERSE
for level = Lvl:-1:1
    for kk = 1:length(hssTL{level})

        Ntmp = hssTL{level}(kk);       % the current Node
        rk   = size( DL{Ntmp},1 );     % the current Node residual rank

        if hssTT{Ntmp}.kids(1) == 0    % base blocks where the data kicks in

            xt{Ntmp} = x( rg(1,Ntmp):rg(2,Ntmp) );   % data kicks in
            if Ntmp == SN
                break;
            end
            xt{Ntmp} = Q1{Ntmp}' * xt{Ntmp};         % Q factor on the right hand side
            nflops = nflops + flops('mv', Q1{Ntmp},'y');

            if ~isempty( L{Ntmp} )  % excluding full rank case: no compression at all
                
                xt{Ntmp}( 1:end-rk ) = L{Ntmp} \ xt{Ntmp} ( 1:end-rk );   % local solution
                nflops = nflops + numel(L{Ntmp});
                
                xt{Ntmp}( end-rk+1:end ) = xt{Ntmp}( end-rk+1:end ) - DL{Ntmp} * xt{Ntmp}( 1:end-rk );   % local update
                nflops = nflops + flops('mv', DL{Ntmp},'n');
                nflops = nflops + flops('sum', xt{Ntmp}( end-rk+1:end ));
                
                Vx{Ntmp} = V{Ntmp}( 1:end-rk,: )' * xt{Ntmp}( 1:end-rk );   % similar to two vectors multiplication
                nflops = nflops + flops('mv', V{Ntmp}( 1:end-rk,: ),'y');
                
            end

        else

            kid1 = hssTT{Ntmp}.kids(1);
            kid2 = hssTT{Ntmp}.kids(2);
            rk1  = size( DL{kid1},1 );
            
            if kid2 == SN 
                
                if ~isempty( Vx{kid1} )
                    xt{kid2} = xt{kid2} - U{kid2} * ( B{kid2} * Vx{kid1} );
                    nflops = nflops + flops('mv',B{kid2},'n');
                    nflops = nflops + flops('mv',U{kid2},'n');
                    nflops = nflops + flops('sum', xt{kid2});
                end
                
                xt{Ntmp} = xt{kid1}( end-rk1+1:end );
                
                if ~isempty( L{Ntmp} )
                     
                    xt{Ntmp} = U{Ntmp} \ (L{Ntmp} \ (DL{Ntmp} * xt{Ntmp}));   % Here L is a square matrix
                    nflops = nflops + numel(L{Ntmp});
                    
                end
                
                xt{kid2} = xt{kid2} - U{kid2} * ( B{kid2} * ( V{kid1}( end-rk1+1:end,: )' * xt{Ntmp} ) );
                nflops = nflops + flops('mv',U{kid2},'n') + flops('mv',B{kid2},'n') + flops('mv',V{kid1}( end-rk1+1:end,: ),'y');
                nflops = nflops + flops('sum', xt{kid2});
                
                break;
                
            end
            
            rk2  = size( DL{kid2},1 );

            % local update within the current Node
            if ~isempty( Vx{kid2} )
                
                xt{kid1}( end-rk1+1:end ) = xt{kid1}( end-rk1+1:end ) - U{kid1} * ( B{kid1} * Vx{kid2} );
                nflops = nflops + flops('mv',B{kid1},'n');
                nflops = nflops + flops('mv',U{kid1},'n');
                nflops = nflops + flops('sum', xt{kid1}( end-rk1+1:end ));
                
            end

            if ~isempty( Vx{kid1} )
                
                xt{kid2}( end-rk2+1:end ) = xt{kid2}( end-rk2+1:end ) - U{kid2} * ( B{kid2} * Vx{kid1} );
                nflops = nflops + flops('mv',B{kid2},'n');
                nflops = nflops + flops('mv',U{kid2},'n');
                nflops = nflops + flops('sum', xt{kid2}( end-rk2+1:end ));
                
            end

            xt{Ntmp} = [ xt{kid1}( end-rk1+1:end ); xt{kid2}( end-rk2+1:end ) ];
            xt{Ntmp} = Q1{Ntmp}' * xt{Ntmp};
            nflops = nflops + flops('mv',Q1{Ntmp},'y');

            if ~isempty( L{Ntmp} )

                xt{Ntmp}( 1:end-rk ) = L{Ntmp} \ xt{Ntmp} ( 1:end-rk );
                nflops = nflops + numel(L{Ntmp});
                
                xt{Ntmp}( end-rk+1:end ) = xt{Ntmp}( end-rk+1:end ) - DL{Ntmp} * xt{Ntmp}( 1:end-rk );
                nflops = nflops + flops('mv',DL{Ntmp},'n');
                nflops = nflops + flops('sum', xt{Ntmp}( end-rk+1:end ));
                
                %%%%%  MAGIC FORMULA BY Jianlin: to accumulate corrections for RHS  %%%%%
                Vx{Ntmp} = V{Ntmp}( 1:end-rk,: )' * xt{Ntmp}( 1:end-rk );
                nflops = nflops + flops('mv',V{Ntmp}( 1:end-rk,: ),'y');

                if ~isempty( Vx{kid1} ) 
                    Vx{Ntmp} = Vx{Ntmp} + W{kid1}' * Vx{kid1};
                    nflops = nflops + flops('mv',W{kid1},'y');
                    nflops = nflops + flops('sum', Vx{Ntmp});
                end
                
                if ~isempty( Vx{kid2} )
                    Vx{Ntmp} = Vx{Ntmp} + W{kid2}' * Vx{kid2};
                    nflops = nflops + flops('mv',W{kid2},'y');
                    nflops = nflops + flops('sum', Vx{Ntmp});
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

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