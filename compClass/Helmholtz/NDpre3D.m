function  [ TT, Tlabel, label, NB, hssTT, hssTL, rg ] = NDpre3D( dim, tos, Lvl, Lswitch, Nx, Ny, Nz )


% 1. elimination tree construction based on the total number of nodes
TT = TT_generate(2^(Lvl+1)-1);


% 2. TL built up: arrage nodes by level
TL = TL_generate(TT);  


% 3. pts: physical span of each subdomain
pts = pts_generate3D( tos, Nx, Ny, Nz, TT, TL ); clear TL;


% 4. Tlabel: span of each subdomain in the global matrix
Tlabel = Tlabel_generate3D( dim, pts );


% 5. SubTree: the subtree of each 2D separator
ST = ST_generate3D(TT, pts); 


% 6. label: reordering of the mesh
label = label_generate3D( Nx, Ny, Nz, pts, TT, ST );


% 7. NB: all neighboring information
if 1 > 0 
    NB = NB_generate3D_label( dim, tos, TT, pts, label );
else
    NB = NB_generate3D( dim, tos, TT, pts, Tlabel, ST );
end
clear pts ST;


% 8. mt
mt = mt_generate3D( TT, Tlabel, NB );


% 9. hssTT, hssTL, rg: HSS structure
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