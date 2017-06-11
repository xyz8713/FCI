function  [ TT, Tlabel, NB, label, hssTT, hssTL, rg ] = NDpre2D( dim, tos, Lvl, Lswitch, Nz, Nx )


% 1. elimination tree construction based on the total number of nodes
TT = TT_generate(2^(Lvl+1)-1);


% 2. TL built up: arrage nodes by level
TL = TL_generate(TT);


% 3. pts: physical span of each subdomain
pts = pts_generate2D( tos, Nz, Nx, TT, TL );  clear TL;


% 4. Tlabel: span of each subdomain in the global matrix
Tlabel = Tlabel_generate2D( dim, pts );


% 5. NB: all neighboring information
NB = NB_generate2D( dim, tos, TT, pts, Tlabel );


% 6. label: reordering of the mesh
label = label_generate2D( Nz, Nx, TT, pts );  clear pts;


% 7. mt
mt = mt_generate2D( TT, Tlabel, NB );


% 8. hssTT, hssTL, rg: HSS structure
hssTT = cell(1,length(TT));
hssTL = cell(1,length(TT));
rg = cell(1,length(TT));

for Ntmp = 1:length(TT)-1
    if TT{Ntmp}.LL <= Lswitch
        [ hssTT{Ntmp}, hssTL{Ntmp}, rg{Ntmp} ] = mat2hssParpre( mt{Ntmp}, 1 );
    end
end
if Lswitch ~= 0
    [ hssTT{length(TT)}, hssTL{length(TT)}, rg{length(TT)} ] = mat2hssParpre( mt{length(TT)}, 0 );
end


return;
end