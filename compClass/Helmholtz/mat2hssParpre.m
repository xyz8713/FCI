function [ hssTT, hssTL, rg ] = mat2hssParpre( m, SchurInd )

% This is the preprocessing step for HSS compression, factorization and
% solution

lm = length(m);
hssTT = TT_generate( 2*lm-1 );

% with Schur to modify hssTT
if SchurInd ~= 0
    
    TT = TT_generate( 2*lm-3 );
    for kk = 1:2*lm-3
        hssTT{kk} = TT{kk};
        hssTT{kk}.LL = hssTT{kk}.LL + 1;
    end
    clear TT;
    hssTT{ 2*lm-3 }.parent = 2*lm-1;
    hssTT{ 2*lm-2 }.parent = 2*lm-1;
    hssTT{ 2*lm-2 }.LL = 2;
    hssTT{ 2*lm-2 }.kids = zeros(1,2);
    hssTT{ 2*lm-1 }.kids(1) = 2*lm-3;
    hssTT{ 2*lm-1 }.kids(2) = 2*lm-2;
  
end

hssTL = TL_generate( hssTT );
rg    = rg_generate( m, hssTT );

return;
end