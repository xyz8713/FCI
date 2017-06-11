%%addpath /project/scicom/scicom00/DATA/MATRICES/Bomhof

addpath /project/scicom/scicom00/DATA/MATRICES/MATLAB_HB

if (1)
    %%  load('psmigr_2');
  load('sherman5') 
    B = Problem.A;
    n = size(B,1);
    if (1)
        d = max(abs(B)) ;
        for j=1:n
            if (d(j) <0), d(j) = -d(j);,end
            B(:,j) = B(:,j)/d(j);
        end
    end
    ex = ones(n,1);
    rhs = B*ex;
else
    load('circuit_3');
    B = Problem.A;
    rhs = Problem.b;
    ex = Problem.x;
    n = size(B,1);
    %    d = ones(n,1);
    %%-------------------- 
    d = max(abs(B)) ;
    for j=1:n
        if (d(j) <0), d(j) = -d(j);,end
        B(:,j) = B(:,j)/d(j);
    end
end