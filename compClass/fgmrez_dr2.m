function [sol,res,its] = fgmrez_dr2 (A,precClass,precfun,rhs,sol,ITopts)
%-----------complex version-----------------------------------------
% function [sol,res,its] = fgmres (A,PRE,rhs,sol)
% restarted gmres with Krylov subspace of dim = im.
% NOTE: this is actually fllexibe (FGMRES) -- allows
% variations in preconditioner
% The first inner FGMRES method, return a initial deflation
% subspace u
%------------------------------------------------------------------
%%

tolIts = ITopts.tolIts;
maxits = ITopts.maxits;
outputG = ITopts.outputG;
im  = ITopts.im;
nc = precClass.nC ;
Ndef = 0;
if (ITopts.Nvec) 
    U = precClass.U;
    Ndef= size(U,2);
end
%fprintf(1,'Ndef = %d\n',Ndef);
%%%% im1 is im+10 injected eigenvectors
im1 = im+Ndef;
%%maxits = im1;

n = size(A,1)    ;
its = 0 ;
%
% main loop
%
sigma = ITopts.shift;

while (its < maxits)
    vv(1:n,1) = rhs - A*sol + sigma*sol  ;
    %%
    ro = norm(vv(1:n,1),2)  ;
    res(its+1) = ro;
    if (its  == 0)
        tol1=tolIts*ro  ;
    end  ;
    if (ro <= tol1 | its >= maxits)
        return
    end
    t = 1.0/ ro;
    vv(1:n,1) = vv(1:n,1) * t  ;
    %       initialize 1-st term  of rhs of hessenberg system..  ;
    rs(1) = ro  ;
    %%-------------------- print its/residual info
    if (outputG)
        fprintf(1,' its %d  res %e \n',its,ro)
    end
    i = 0  ;
    %-------------------- inner gmres loop
    while (i < im1  &  (ro  >  tol1)  &  its < maxits)
        i=i+1  ;
        its = its + 1  ;
        i1 = i + 1 ;
        %%-------------------- for first vectors Use vv vector
        if (i<=Ndef)
            w = U(:,i);
        else
        %%-------------------- else use deflation vector
            w = vv(:,i) ; 
        end

        z = feval(precfun,precClass,w);
        W(:,i) = z;
        
        z = A*z - sigma*z;
        %%z = A*z ;
        %%--------------------       modified GS  ;
        for j=1:i
            t = vv(1:n,j)'*z  ;
            hh(j,i) = t  ;
            z = z - t*vv(1:n,j)  ;
        end  
        t = norm(z,2)  ;
        hh(i1,i) = t  ;
        if (t  ~= 0.0)
            t = 1.0 / t  ;
            vv(1:n,i1) = z*t  ;
        end   %% IF
        %%
        if (i ~= 1)
            %
            %--------------------previous rots. on i-th column of h  ;
            %
            for k=2:i
                k1 = k-1  ;
                t = hh(k1,i)  ;
                hh(k1,i) = conj(c(k1))*t + s(k1)*hh(k,i)  ;
                hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)  ;
            end  %% FOR
        end   %% IF
        %%
        gam = sqrt(abs(hh(i,i))^2 + abs(hh(i1,i))^2)  ;
        if (gam  == 0.0)
            gam = eps;
        end  
        %
        %       determine plane rotation and update rhs of ls pb   ;
        %
        c(i) = hh(i,i)/gam  ;
        s(i) = hh(i1,i)/gam  ;
        rs(i1) = -s(i)*rs(i)  ;
        rs(i) =  conj(c(i))*rs(i)  ;
        %-------------------- test for convergence-  ;
        hh(i,i) = conj(c(i))*hh(i,i) + s(i)*hh(i1,i)  ;
        ro = abs(rs(i1))  ;
        sprintf(' ro = \n',ro) ;
        res(its+1) = ro   ;
        if (outputG)
            fprintf(1,' its %d  res %e \n',its,ro)
        end
    end    %% end of while (im) loop
    %
    %--------------------now compute solution. 
    %
    rs(i) = rs(i)/hh(i,i)  ;
    for  k=i-1:-1:1  ;
         t=rs(k)  ;
         for j=k+1:i  ;
            t = t-hh(k,j)*rs(j)  ;
         end  ;
         rs(k) = t/hh(k,k)  ;
    end  
    %    done with back substitution..  ;
    %    now form linear combination to get solution  ;
    %       
    rs = rs(:);
    temp = W(:,1:i)*rs(1:i);
    temp =  feval(precfun,precClass,temp);
    sol = sol + temp;
    if ((ro  <=  tol1) | (its >= maxits))
        fprintf(1,' total its %d  final rel res %e \n',its,ro/abs(res(1)))
        if (ITopts.Nvec == 0) 
            return;
        end
        %%-------------------- ritz values version
        HH = hh(1:i,1:i);
        GG = vv(:,1:i)'*W(:,1:i);
        [Vec,Val] = eig(HH,GG);
        [~,I] = sort(abs(diag(Val)-sigma), 'ascend' );
        Ndef = min(ITopts.Nvec,i);
        Q  = Vec(:,I(1:Ndef));
        U = W(:,1:i)*Q;
        %%-------------------- end
        precClass.U = U;
        return;
    end
end  %% end while -- restart
end


