function  [ x, nflops ] = mfsolvePar( x, Lswitch, TT, Tlabel, NB, hssTT, hssTL, rg, Q1, Q2, L, DL, U, B, W, V )

% the entire multifrontal solver is also parallel. We traverse the entire
% elimination tree by different levels

% flops of solver calculation
nflops = 0;

%%%%%%%%%%%%%%%%%   L y = b   %%%%%%%%%%%%%%%%%%%%%
% from bottom to top

for Ntmp = 1:length(TT)
        
    Nhead = Tlabel(Ntmp);         % current Node head
    Ntail = Tlabel(Ntmp+1) - 1;   % current Node tail

    % full LU solver, no HSS structure for L > Lswitch
    if TT{Ntmp}.LL > Lswitch

        x( Nhead:Ntail ) = L{Ntmp} \ x( Nhead:Ntail );
        nflops = nflops + numel(L{Ntmp});

        counter = 1;
        for ii = 1:length(NB{Ntmp})
            Phead = NB{Ntmp}{ii}.Ohead;
            Ptail = NB{Ntmp}{ii}.Otail;
            x( Phead:Ptail ) = x( Phead:Ptail ) - DL{Ntmp}( counter:counter+Ptail-Phead,:) * x( Nhead:Ntail );
            nflops = nflops + flops('mv', DL{Ntmp}( counter:counter+Ptail-Phead, :),'n');
            nflops = nflops + flops('sum', x( Phead: Ptail ));
            counter = counter + Ptail - Phead + 1;
        end

    % HSS solver
    else

        xt = zeros( rg{Ntmp}(2,end),1 );   % the size of the current frontal matrix
        xt( 1:Ntail-Nhead+1 ) = x( Nhead:Ntail );

        if Ntmp ~= length(TT)

            % copy Schur complement RHS
            counter = Ntail - Nhead + 2;
            for ii = 1:length(NB{Ntmp})
                Phead = NB{Ntmp}{ii}.Ohead;
                Ptail = NB{Ntmp}{ii}.Otail;
                xt( counter:counter+Ptail-Phead ) = x( Phead:Ptail );
                counter = counter + Ptail - Phead + 1;
            end

            % partial solver for HSS
            [ xt, nflops1 ] = mat2hssParsolL( xt,Q1{Ntmp},L{Ntmp},DL{Ntmp},U{Ntmp},B{Ntmp},W{Ntmp},V{Ntmp},hssTT{Ntmp},hssTL{Ntmp},rg{Ntmp} );
            nflops = nflops + nflops1; clear nflops1;

        else

            [ xt, nflops1 ] = mat2hssParsol( xt,Q1{Ntmp},Q2{Ntmp},L{Ntmp},DL{Ntmp},U{Ntmp},B{Ntmp},W{Ntmp},V{Ntmp},hssTT{Ntmp},hssTL{Ntmp},rg{Ntmp} );
            nflops = nflops + nflops1; clear nflops1;

        end

        % copy back to the original solution x
        x( Nhead:Ntail ) = xt( 1:Ntail-Nhead+1 );
        if Ntmp ~= length(TT)
            counter = Ntail - Nhead + 2;
            for ii = 1:length(NB{Ntmp})
                Phead = NB{Ntmp}{ii}.Ohead;
                Ptail = NB{Ntmp}{ii}.Otail;
                x( Phead:Ptail ) = xt( counter:counter+Ptail-Phead );
                counter = counter + Ptail - Phead + 1;
            end
        end
        clear xt; 

    end
    
end


%%%%%%%%%%%%%%%%%   U x = y   %%%%%%%%%%%%%%%%%%%%%
% from top to bottom

for Ntmp = length(TT):-1:1
           
    Nhead = Tlabel(Ntmp);         % current Node head
    Ntail = Tlabel(Ntmp+1) - 1;   % current Node tail

    if TT{Ntmp}.LL > Lswitch
        
        if Ntmp ~= length(TT)
            counter = 1;
            for ii = 1:length(NB{Ntmp})
                Phead = NB{Ntmp}{ii}.Ohead;
                Ptail = NB{Ntmp}{ii}.Otail;
                x( Nhead:Ntail ) = x( Nhead:Ntail ) - V{Ntmp}( :, counter:counter+Ptail-Phead ) * x( Phead:Ptail );
                nflops = nflops + flops('mv', V{Ntmp}( :, counter:counter+Ptail-Phead ),'n');
                nflops = nflops + flops('sum', x( Nhead:Ntail ));
                counter = counter + Ptail - Phead + 1;
            end
        end
        x( Nhead:Ntail ) = U{Ntmp} \ x( Nhead:Ntail );
        nflops = nflops + numel(U{Ntmp});

    else

        if Ntmp ~= length(TT)
            
            xt = zeros( rg{Ntmp}(2,end),1 );
            xt( 1:Ntail-Nhead+1 ) = x( Nhead:Ntail );

            counter = Ntail - Nhead + 2;
            for ii = 1:length(NB{Ntmp})
                Phead = NB{Ntmp}{ii}.Ohead;
                Ptail = NB{Ntmp}{ii}.Otail;
                xt( counter:counter+Ptail-Phead ) = x( Phead:Ptail );
                counter = counter + Ptail - Phead + 1;
            end

            [ xt, nflops1 ] = mat2hssParsolU( xt,Q2{Ntmp},L{Ntmp},DL{Ntmp},U{Ntmp},B{Ntmp},V{Ntmp},hssTT{Ntmp},hssTL{Ntmp},rg{Ntmp} );
            nflops = nflops + nflops1; clear nflops1;
            x( Nhead:Ntail ) = xt( 1:Ntail-Nhead+1 );
            clear xt;
        
        end

    end
    
end

return;
end