function [P, E, rk] = SRRQR_Wang( A, F, tol )
% strong rank revealing QR using Ming Gu's paper in 1996
% U = [eye(rk);E];
% U(P,:) = U;
% A = U * A(P(1:rk),:);
% F here is the ratio of the determinant, 1.5 < F < 2
% tol is for RRQR

% conventional RRQR of A.', and estimating the rank
[Q,R,P] = qr(A.',0);
N = size(A,1);
M = sum(abs(diag(R)) >= abs(R(1,1)) * tol);    % M is the rank
R = R(1:M,:); clear Q;

% to judge M and N
if M == N
    E = [];
else
    while 1 > 0
        
        % find the current rho: max(abs( A^{-1}B ))
        E = R(:,1:M) \ R(:,M+1:N);
        [vm,im] = max(abs(E));    % caution: abs value
        [rho,K] = max(vm);
        I       = im(K);
        K       = K + M;
        
        % to judge
        % display(['rho = ',num2str(rho)]);
        if rho <= F
            break;
        else
            % interchange the column M and K
            R(:,[M,K]) = R(:,[K,M]);
            P(  [M,K]) = P(  [K,M]);
            % I == M no need to do Givens transformation
            if I < M
                % interchange the column I and M
                R(:,[I,M]) = R(:,[M,I]);
                P(  [I,M]) = P(  [M,I]); 
                % upper Hessenburg matrix
                for ii = M:-1:I+2    % caution: I+2
                    x1 = R(ii-1,I);
                    x2 = R(ii  ,I);
                    R(ii-1,I) = sqrt(x1^2+x2^2);
                    R(ii,  I) = 0;
                    if R(ii-1,I) ~= 0
                        G = [x1,x2;-x2,x1]/R(ii-1,I);
                        R(ii-1:ii,ii-1:N) = G * R(ii-1:ii,ii-1:N);
                    end
                end
                % upper triangular matrix
                for ii = I:M-1
                    x1 = R(ii  ,ii);
                    x2 = R(ii+1,ii);
                    R(ii  ,ii) = sqrt(x1^2+x2^2);
                    R(ii+1,ii) = 0;
                    if R(ii,ii) ~= 0
                        G = [x1,x2;-x2,x1]/R(ii,ii);
                        R(ii:ii+1,ii+1:N) = G * R(ii:ii+1,ii+1:N);
                    end
                end
            end
        end
        
    end
end

E = E.';  rk = M;

return;
end