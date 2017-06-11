function [ Q1, Q2, L, DL, U, B, W, V, S, nflops, storage ] = mat2hssParfac( D, U, R, B, W, V, hssTT, hssTL, SchurInd )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  HSS parallel factorization  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if SchurInd == 0
%     S = [];
% else
%     S = Schur complement;
% end

% nflops and storage initialization
nflops  = 0;
storage = 0;

lenTT = length(hssTT);   % length of the tree
Lvl   = length(hssTL);   % total level

Q1 = cell(1,lenTT);   % row multiplier orthogonal matrix
Q2 = cell(1,lenTT);   % col multiplier orthogonal matrix
L  = cell(1,lenTT);   % diagonal lower triangular matrix
DL = cell(1,lenTT);   % diagonal residual matrix

%%%%%%%  NO Schur complement  %%%%%%%
if SchurInd == 0  
    
    S = [];
    %%%%%%%%%%%%  Factorization loop  %%%%%%%%%%%%%%
    for level = Lvl:-1:1
        for kk = 1:length(hssTL{level})

            Ntmp = hssTL{level}(kk);    % the current Node

            if hssTT{Ntmp}.kids(1) == 0

                Di = D{Ntmp};
                Ui = U{Ntmp};
                Vi = V{Ntmp};

            else

                kid1 = hssTT{Ntmp}.kids(1);
                kid2 = hssTT{Ntmp}.kids(2);
                rk1 = size( D{kid1},1 );
                rk2 = size( D{kid2},1 );

                Di = [ D{kid1}, U{kid1} * B{kid1} * V{kid2}( end-rk2+1:end, : )' ;  ...
                       U{kid2} * B{kid2} * V{kid1}( end-rk1+1:end, : )', D{kid2} ] ;
                % m*n, n*l, l*k three matrices multiplication
                % flops = m*((2n-1)*l+(2l-1)*k)
                nflops = nflops + size(U{kid1},1) * ( (2*size(U{kid1},2)-1)*size(B{kid1},2) + (2*size(B{kid1},2)-1)*rk2 );
                nflops = nflops + size(U{kid2},1) * ( (2*size(U{kid2},2)-1)*size(B{kid2},2) + (2*size(B{kid2},2)-1)*rk1 );

                if level > 1
                    Ui = [ U{kid1} * R{kid1}; U{kid2} * R{kid2} ];
                    Vi = [ V{kid1}( end-rk1+1:end, : ) * W{kid1}; V{kid2}( end-rk2+1:end, : ) * W{kid2} ];
                    nflops = nflops + flops('prod',U{kid1},'n',R{kid1},'n') + flops('prod',U{kid2},'n',R{kid2},'n');
                    nflops = nflops + rk1*(2*size(W{kid1},1)-1)*size(W{kid1},2) + rk2*(2*size(W{kid2},1)-1)*size(W{kid2},2);
                end
                
                clear rk1 rk2;

            end
            
            if level > 1

                % the rank of compression: size( Ui, 2 )
                sd = size( Ui,1 );   % size of the current D
                rk = size( Ui,2 );   % rank 

                % 1. QL factorization for Ui
                storage = storage + (sd-(rk-1)/2)*rk;  % storage for Q1 factor
                [ Q1{Ntmp}, U{Ntmp} ] = ql( Ui );
                U{Ntmp} = U{Ntmp}( sd-rk+1:sd, : );
                nflops = nflops + flops('qr', Ui, 'n');

                % 2. update Di
                Di = Q1{Ntmp}' * Di;
                nflops = nflops + sd*rk*(2*sd-rk);

                % 3. QR factorization of the upper block of Di
                storage = storage + (sd-(sd-rk-1)/2)*(sd-rk);   % storage for Q2 factor
                [ Q2{Ntmp}, tmp ] = qr( Di( 1:sd-rk, : )' );
                L{Ntmp} = tmp( 1:sd-rk, : )';  clear tmp;
                nflops = nflops + flops( 'qr', Di( 1:sd-rk, : ),'y' );

                % 4. update V{Ntmp}
                V{Ntmp} = Q2{Ntmp}' * Vi;
                nflops = nflops + size(Vi,2)*(sd-rk)*(sd+rk);

                % 5. update D{Ntmp}, L{Ntmp}, DL{Ntmp}
                nflops = nflops + rk*sd*(2*sd-rk);
                Di = Di( sd-rk+1:sd, : ) * Q2{Ntmp};
                DL{Ntmp} = Di( :, 1:sd-rk );
                D{Ntmp}  = Di( :, sd-rk+1:sd );

            else

                [ L{Ntmp}, U{Ntmp}, DL{Ntmp} ] = lu(Di);    % factorize the final block

            end

            clear Di Ui Vi;
            storage = storage + numel(L{Ntmp}) + numel(DL{Ntmp}) + numel(U{Ntmp}) + numel(B{Ntmp}) + numel(W{Ntmp}) + numel(V{Ntmp});

        end
    end
    
    
else   %%%%%  with Schur complement  %%%%%
    
    SN = lenTT - 1;    % Schur complement Node
    
    %%%%%%%%%%%%  Factorization loop  %%%%%%%%%%%%%%
    for level = Lvl:-1:1
        for kk = 1:length(hssTL{level})

            Ntmp = hssTL{level}(kk);    % the current Node
            
            if Ntmp == SN   % bypass the Schur Node
                break;
            end

            if hssTT{Ntmp}.kids(1) == 0

                Di = D{Ntmp};
                Ui = U{Ntmp};
                Vi = V{Ntmp};

            else

                kid1 = hssTT{Ntmp}.kids(1);
                kid2 = hssTT{Ntmp}.kids(2);
                rk1 = size( D{kid1},1 );
                rk2 = size( D{kid2},1 );

                if kid2 == SN   % which means the root Node 
                    
                    Di = D{kid1};
                    Ui = U{kid1};
                    Vi = V{kid1}( end-rk1+1:end, : );
                    [ L{Ntmp}, U{Ntmp}, DL{Ntmp} ] = lu(Di);

                    tmp = Vi' * (U{Ntmp} \ (L{Ntmp} \ (DL{Ntmp} * Ui)));   % factorize the final block
                    nflops = nflops + 2/3*size(Di,1)^3*size(Ui,2) + flops('prod',Vi,'y',Ui,'n');
                    
                    tmp = B{SN} * tmp * B{SN-1};
                    nflops = nflops + size(B{SN},1) * ( (2*size(B{SN},2)-1)*size(tmp,2) + (2*size(tmp,2)-1)*size(B{SN-1},2) );
                    
                    S = D{SN} - U{SN} * tmp * V{SN}' ;
                    nflops = nflops + size(U{SN},1) * ( (2*size(U{SN},2)-1)*size(tmp,2) + (2*size(tmp,2)-1)*size(V{SN},1) );
                    nflops = nflops + flops('sum',D{SN});
                    clear tmp;
                    
                    break;
                    
                else
                    
                    Di = [ D{kid1}, U{kid1} * B{kid1} * V{kid2}( end-rk2+1:end, : )' ;  ...
                           U{kid2} * B{kid2} * V{kid1}( end-rk1+1:end, : )', D{kid2} ] ;
                    nflops = nflops + size(U{kid1},1) * ( (2*size(U{kid1},2)-1)*size(B{kid1},2) + (2*size(B{kid1},2)-1)*rk2 );
                    nflops = nflops + size(U{kid2},1) * ( (2*size(U{kid2},2)-1)*size(B{kid2},2) + (2*size(B{kid2},2)-1)*rk1 );
                       
                    Ui = [ U{kid1} * R{kid1}; U{kid2} * R{kid2} ];
                    Vi = [ V{kid1}( end-rk1+1:end, : ) * W{kid1}; V{kid2}( end-rk2+1:end, : ) * W{kid2} ];
                    nflops = nflops + flops('prod',U{kid1},'n',R{kid1},'n') + flops('prod',U{kid2},'n',R{kid2},'n');
                    nflops = nflops + rk1*(2*size(W{kid1},1)-1)*size(W{kid1},2) + rk2*(2*size(W{kid2},1)-1)*size(W{kid2},2);
                    
                end
                
                clear rk1 rk2;

            end

            % the rank of compression: size( Ui, 2 )
            sd = size( Ui,1 );   % size of the current D
            rk = size( Ui,2 );   % rank 

            % 1. QL factorization for Ui
            storage = storage + (sd-(rk-1)/2)*rk;  % storage for Q1 factor
            [ Q1{Ntmp}, U{Ntmp} ] = ql( Ui );
            U{Ntmp} = U{Ntmp}( sd-rk+1:sd, : );
            nflops = nflops + flops('qr', Ui,'n');

            % 2. update Di
            Di = Q1{Ntmp}' * Di;
            nflops = nflops + sd*rk*(2*sd-rk);

            % 3. QR factorization of the upper block of Di
            storage = storage + (sd-(sd-rk-1)/2)*(sd-rk);   % storage for Q2 factor
            [ Q2{Ntmp}, tmp ] = qr( Di( 1:sd-rk, : )' );
            L{Ntmp} = tmp( 1:sd-rk, : )';  clear tmp;
            nflops = nflops + flops( 'qr', Di( 1:sd-rk, : ),'y' );

            % 4. update V{Ntmp}
            V{Ntmp} = Q2{Ntmp}' * Vi;
            nflops = nflops + size(Vi,2)*(sd-rk)*(sd+rk);

            % 5. update D{Ntmp}, L{Ntmp}, DL{Ntmp}
            nflops = nflops + rk*sd*(2*sd-rk);
            Di = Di( sd-rk+1:sd, : ) * Q2{Ntmp}; 
            DL{Ntmp} = Di( :, 1:sd-rk );
            D{Ntmp}  = Di( :, sd-rk+1:sd );
            
            clear Di Ui Vi;
            storage = storage + numel(L{Ntmp}) + numel(DL{Ntmp}) + numel(U{Ntmp}) + numel(B{Ntmp}) + numel(W{Ntmp}) + numel(V{Ntmp});

        end
    end
    
end


return;
end