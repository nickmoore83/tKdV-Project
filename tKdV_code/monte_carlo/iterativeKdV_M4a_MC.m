function [zk,ik]=iterativeKdV_M4a_MC(zk0, tau,yk,p, tol)
% Function iterates to solve the solution for Z-Y=tau/2*F(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0 = p.D0; E0 = p.E0; Ls=p.Ls;
J = p.J; MC = p.MC;
alp1 = p.alpha1; alp2 = p.alpha2;

k = repmat([0:J/2 -J/2+1:-1]',[1 MC]);
Lk1 = (2*Ls)^(-3)*D0^(alp2)*1i*k.^3;
M = 1./(1-tau/2*Lk1);
Dk1 = E0^(1/2)*(2*Ls)^(-3/2)*D0^(alp1)*1i*k;
% tol=1e-10; % iteration tolerance

zk1=zk0;
zk2=zk0;
err0=0;
Ni=500;  % maximum interation steps
for ll=1:Ni
    
    zk2=M.*(yk-tau/2*.5*Dk1.*u2k_dealiasing_MC(zk1,p));
    
    err=sum(abs(zk2-zk1).^2);
    if max(err)<tol && sum(err)>=sum(err0)
%         display(['iteration steps N = ',num2str(ll)]);
        ik=0;
        break;
    end
    err0=err;
    zk1=zk2;
    if ll==Ni
        disp(['maximum interation N = ',num2str(ll),' reached!']);
        disp(['err = ',num2str(max(err)),', sum = ', num2str(sum(err>=err0))]);
        err
        ik=1;
    end
end
zk=zk2;