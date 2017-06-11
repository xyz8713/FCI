function [ Q, R, nflops, rk ] = compr_wang( A, tol, ind )

% compression algorithm by Modified Gram-Schimdt Column Pivoting (mgscolpv)
% algorithm, rank-revealing strategy

nflops = 0;

if ind == 0    % using matlab built-in function to carry out QR with pivoting
    
    [Q,R,P] = qr(A);
    rk = sum( abs(diag(R)) >= abs(R(1,1))*tol );
    Q = Q(:,1:rk);
    R = R(1:rk,:);
    R = R*P';
    
else    % self-builtup function

    [m,n] = size(A);
    LDN = min( m,n );   % leading dimension of rank 

    R   = zeros( LDN,n );
    piv = zeros( 1,LDN );    % pivoting vector: Lapack rule

    % initialize col_norm: square of norms of each column
    col_norm = zeros(1,n);
    for kk = 1:n
        col_norm(kk) = norm(A(:,kk));   % IMPORTANT: norm
    end
    nflops = nflops + n*(2*m-1);

    for kk = 1:LDN

        % searching for the largest value of norms 
        [ vmax, piv(kk) ] = max( col_norm(kk:n) );
        piv(kk) = piv(kk) + kk-1 ;

        % colum pivoting A and R: IMPORTANT!!!
        if piv(kk) > kk
            A(:,[kk,piv(kk)]) = A(:,[piv(kk),kk]);
            R(1:kk-1,[kk,piv(kk)]) = R(1:kk-1,[piv(kk),kk]);
            col_norm(piv(kk)) = col_norm(kk);
        end

        % judging whether to continue or not
        R(kk,kk) = vmax;
        if abs(R(kk,kk)) < tol * abs(R(1,1))
            rk = kk-1;
            rk_piv = kk;
            break;
        else
            rk = kk;
            rk_piv = kk;
        end

        % orthonormalization and updating
        A(:,kk) = A(:,kk) / R(kk,kk);
        nflops = nflops + m;

        R(kk,kk+1:n) = A(:,kk)' * A(:,kk+1:n);
        nflops = nflops + (n-kk)*(2*m-1);

        A(:,kk+1:n) = A(:,kk+1:n) - A(:,kk) * R(kk,kk+1:n);
        nflops = nflops + (n-kk)*2*m;

        for ii = kk+1:n
            col_norm(ii) = norm(A(:,ii));  % CAREFUL: norm instead of update
        end
        nflops = nflops + (n-kk)*(2*m-1);

    end

    piv = piv(1:rk_piv);
    Q = A(:,1:rk);
    R = R(1:rk,:); 
    % applying the transpose of permutation: A P = Q * R, A = Q * R * P'
    for kk = rk_piv:-1:1
        if piv(kk) > kk
            R(:,[kk,piv(kk)]) = R(:,[piv(kk),kk]);
        end
    end

end

return;
end