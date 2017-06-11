function [ Q, L ] = ql( A )

% A S = Q P P L S

[ Q, L ] = qr( full(A(:,end:-1:1)) );
Q = Q(:,end:-1:1);
L = L(end:-1:1,end:-1:1);

return;
end