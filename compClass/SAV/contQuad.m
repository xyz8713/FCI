    function [zk, omega ] =  contQuad(n , method) 
%%  function [zk, omega ] =  contQuad(n , method) 
%% method: 1= trapz, 2 = mid-pt, 0 - Gauss.

if( method == 0 )
    disp 'Use Gauss Legendre' ;
%-------------------- Gauss-Legendre
    n = n - 1 ; % n + 1 points without this line
    beta = .5 ./ sqrt(1 - (2 * ( 1 : n ) ).^( - 2 ) );
    T = diag( beta, 1 ) + diag( beta, - 1 ) ;
    [ V, D ] = eig( T ) ;
    x = diag( D ); 
    [ x, I ] = sort( x );
    omega = ( V( 1, I ).^2 )';
    theta = (pi / 2) .* (1 - x) ;
elseif( method == 1 )
    disp 'Use Gauss Jacobi';
    [ x, omega ] = gauss_jacobi( n, 0, 0 );
    theta = pi*x;
else
    disp 'midpoint' ;
    theta = pi *(2*[1:n]-1)/(2*n);
    omega = ones(1,n)/n;
end
%%-------------------- 
 zk    = exp(1i * theta);
%%-------------------- In order to use standard form 
%%                     sum_k om_l /(z-sig_k) 
%%                     need to multiply by -1/4
%% a factor of 1/2 comes from the integration formula
%% changing from z variable on circle to x variable in [-1, 1]
%% the other half comes from the fact tha the h should be halved
%% since the number of points is 2*n not n (the other half-circle).
%omega = -omega .* zk/2.00 ; 
 omega = -omega  .* zk/2;
%% 
