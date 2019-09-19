
function [zk,ik]=iterativeKdV_M4a_MC(zk0, tau,yk,p, tol)
% Function iterates to solve the solution for Z-Y=tau/2*F(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2=p.C2; C3=p.C3;
D0=p.D0;
J = p.J; MC = p.MC;

k = repmat([0:J/2 -J/2+1:-1]',[1 MC]);
Lk1 = C2*D0^(1/2)*1i*k.^3;
M = 1./(1-tau/2*Lk1);
Dk1 = C3*D0^(-3/2)*1i*k;
% tol=1e-10; % iteration tolerance

zk1=zk0;
zk2=zk0;
err0=0;
Ni=500;  % maximum interation steps
for ll=1:Ni
    
    zk2=M.*(yk-tau/2*.5*Dk1.*u2k_dealiasing_MC(zk1));
    
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
%         display(err)
        ik=1;
    end
end
zk=zk2;