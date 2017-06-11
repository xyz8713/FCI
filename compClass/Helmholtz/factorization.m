function  [ Q1, Q2, L, DL, U, B, W, V, nflops, storage ] = factorization( A, tol, Lswitch, TT, Tlabel, NB, hssTT, hssTL, rg )

% nested dissection based LU decomposition: using extend and add algeorithm
% in multifrontal method
% for level > Lswitch, Lcorner and Lleft are stored in L{Ntmp} and DL{Ntmp},
% Ucorner and Uright are stored in U{Ntmp} and V{Ntmp}

% factorization flops and storage calculation initialization
nflops  = 0; 
storage = 0;

% Stack structure to store all dense matrices
stackInd = 0;
stack = {};       % stores updated matrices waiting to be extend-add

Q1 = cell(1,length(TT));
Q2 = cell(1,length(TT));
L  = cell(1,length(TT));
DL = cell(1,length(TT));
U  = cell(1,length(TT));
B  = cell(1,length(TT));
W  = cell(1,length(TT));
V  = cell(1,length(TT));


for Ntmp = 1:length(TT)     % outer loop for all Nodes
    

    display(['----- the current Node is ',num2str(Ntmp),', among ',num2str(length(TT)),' Nodes. -----']);

    Nhead = Tlabel(Ntmp);          % current node head
    Ntail = Tlabel(Ntmp+1)-1;      % current node tail
    
    %%%%%%%%%%%%%%%  copy matrix A  %%%%%%%%%%%%%%%%%
    sz1 = Ntail - Nhead + 1;
    L{Ntmp} = full(A(Nhead:Ntail,Nhead:Ntail));
    
    if Ntmp ~= length(TT)
        
        sz2 = 0;
        for ii = 1:length(NB{Ntmp})
            sz2 = sz2 + (NB{Ntmp}{ii}.Otail - NB{Ntmp}{ii}.Ohead + 1);
        end
        DL{Ntmp} = zeros(sz2,sz1);
        V{Ntmp}  = zeros(sz1,sz2);
        Schur    = zeros(sz2,sz2);
        
        counter = 1;
        for ii = 1:length(NB{Ntmp})
            Phead = NB{Ntmp}{ii}.Ohead;
            Ptail = NB{Ntmp}{ii}.Otail;
            DL{Ntmp}(counter:counter+Ptail-Phead,1:sz1) = full(A(Phead:Ptail,Nhead:Ntail));
            V{Ntmp}(1:sz1,counter:counter+Ptail-Phead)  = full(A(Nhead:Ntail,Phead:Ptail));
            counter = counter + Ptail-Phead+1;
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  extend-add part for non-lowest level nodes  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if TT{Ntmp}.kids(1) ~= 0   % except the lowest level: non-leaf nodes
        
        if Ntmp == length(TT)
            
            L{Ntmp} = L{Ntmp} + stack{stackInd} + stack{stackInd-1};
            clear stack;
            
        else
        
            % Top Update matrix poped up from the stack
            [ L{Ntmp}, V{Ntmp}, DL{Ntmp}, Schur, nflops1 ] = extadd( L{Ntmp}, V{Ntmp}, DL{Ntmp}, Schur, NB{Ntmp}, stack{stackInd}, NB{TT{Ntmp}.kids(2)} );
            nflops = nflops + nflops1; clear nflops1;
            stack{stackInd}   = [];     % pushed out of the stack
            stackInd = stackInd - 1;


            % Bot Update matrix poped up from the stack
            [ L{Ntmp}, V{Ntmp}, DL{Ntmp}, Schur, nflops1 ] = extadd( L{Ntmp}, V{Ntmp}, DL{Ntmp}, Schur, NB{Ntmp}, stack{stackInd}, NB{TT{Ntmp}.kids(1)} );
            nflops = nflops + nflops1; clear nflops1;
            stack{stackInd}   = [];     % pushed out of the stack
            stackInd = stackInd - 1;
        
        end
        
    end
    
%     DM = [ L{Ntmp},V{Ntmp} ; DL{Ntmp},Schur ]; DMsize = size(DM,1);
%     spy(DM); title(['Node = ',num2str(Ntmp),' among ',num2str(length(TT)),' nodes']); hold on;
%     plot((Ntail-Nhead+1)*ones(1,Ntail-Nhead+1),1:Ntail-Nhead+1,'r--','LineWidth',3); hold on;
%     plot((Ntail-Nhead+1)*ones(1,DMsize-(Ntail-Nhead+1)),Ntail-Nhead+2:DMsize,'g--','LineWidth',3); hold on;
%     plot(1:Ntail-Nhead+1,(Ntail-Nhead+1)*ones(1,Ntail-Nhead+1),'r--','LineWidth',3); hold on;
%     plot(Ntail-Nhead+2:DMsize,(Ntail-Nhead+1)*ones(1,DMsize-(Ntail-Nhead+1)),'g--','LineWidth',3); hold off;
%     pause;
    
    
    %%%%%%%%%%%%%%%  HSS compression  %%%%%%%%%%%%%%%%%
    
    if TT{Ntmp}.LL <= Lswitch    % HSS compression above switch level
        
        if Ntmp == length(TT)
            SchurInd = 0;
        else
            SchurInd = 1;
        end
        
        stackInd = stackInd + 1;
        DM = [ L{Ntmp},V{Ntmp} ; DL{Ntmp},Schur ];
        
        [ Q1{Ntmp},Q2{Ntmp},L{Ntmp},DL{Ntmp},U{Ntmp},B{Ntmp},W{Ntmp},V{Ntmp},stack{stackInd},nflops1,storage1 ] ... 
        = mat2hssPar( DM, hssTT{Ntmp}, hssTL{Ntmp}, rg{Ntmp}, tol, SchurInd );
        
        nflops = nflops + nflops1; clear nflops1;
        storage = storage + storage1; clear storage1; 
            
        clear DM;
        
    else   % full LU decomposition above switch level
        
        [ L{Ntmp} , U{Ntmp} ] = lu( L{Ntmp} ); 
        nflops = nflops + flops('lu', L{Ntmp});
        storage = storage + (Ntail-Nhead+1)^2;    % we use the style of Lapack to store L and U factors
        
        if Ntmp ~= length(TT)
            DL{Ntmp} = DL{Ntmp} / U{Ntmp};
            V{Ntmp}  = L{Ntmp} \ V{Ntmp};
            stackInd = stackInd + 1;
            stack{stackInd} = Schur - DL{Ntmp} * V{Ntmp};

            nflops = nflops + size(U{Ntmp},1)^2 * sz2;
            nflops = nflops + size(L{Ntmp},1)^2 * sz2;
            nflops = nflops + flops('prod', DL{Ntmp}, 'n', V{Ntmp}, 'n');
            nflops = nflops + flops('sum', Schur);
            storage = storage + numel(DL{Ntmp}) + numel(V{Ntmp});
        end
        
    end

end 


return;
end