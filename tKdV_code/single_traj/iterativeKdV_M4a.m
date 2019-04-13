function [zk,ik]=iterativeKdV_M4a(zk0, tau,yk,p, tol)
% Function iterates to solve the solution for Z-Y=tau/2*F(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D0 = p.D0; E0 = p.E0; Ls=p.Ls;
alp1 = p.alpha1; alp2 = p.alpha2;

k = [0:p.J/2 -p.J/2+1:-1]'; % wavenumbers
Lk = (2*Ls)^(-3)*D0^(alp2)*1i*k.^3;
Dk = E0^(1/2)*(2*Ls)^(-3/2)*D0^(alp1)*1i*k;


zk1=zk0;
zk2=zk0;
err0=0;
Ni=1000;  % maximum interation steps
for ll=1:Ni
    M = 1./(1-tau/2*Lk);
    zk2=M.*(yk-tau/2*.5*Dk.*u2k_dealiasing(zk1,p));
    
    err=sum(abs(zk2-zk1).^2);
    if err<tol && err>=err0
%         display(['iteration steps N = ',num2str(ll)]);
        ik=0;
        break;
    end
    err0=err;
    zk1=zk2;
    if ll==Ni
        disp(['maximum interation N = ',num2str(ll),' reached!']);
        ik=1;
    end
end
zk=zk2;