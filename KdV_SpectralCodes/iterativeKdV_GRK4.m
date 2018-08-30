function [uk1,uk2,ik]=iterativeKdV_GRK4(uk10,uk20, uk,p)
% Function iterates to solve the solution for Z-Y=tau/2*F(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=p.dt; A=p.A;

k = (pi/p.Ls)*[0:p.J/2 -p.J/2+1:-1]'; % wavenumbers
Lk = 1i*k.^3;
Dk = 1i*k;
tol=5e-10; % iteration tolerance

uk1=uk10; uk1p=uk1;
uk2=uk20; uk2p=uk2;
err1p=0; err2p=0;
Ni=100;  % maximum interation steps
for ll=1:Ni
    u2k1 = u2k_dealiasing(uk1,p);
    u2k2 = u2k_dealiasing(uk2,p);
    M1 = 1./((1-dt*A(1,1).*Lk).*(1-dt*A(2,2).*Lk) - (dt*A(2,1).*Lk).*(dt*A(1,2).*Lk));
    M2 = 1./((1-dt*A(2,2).*Lk).*(1-dt*A(1,1).*Lk) - (dt*A(1,2).*Lk).*(dt*A(2,1).*Lk));
    
    uk1 = M1.*(uk-3*dt*A(1,1).*Dk.*u2k1-3*dt*A(1,2).*Dk.*u2k2).*(1-dt*A(2,2).*Lk) ...
         +M1.*(uk-3*dt*A(2,1).*Dk.*u2k1-3*dt*A(2,2).*Dk.*u2k2).*(dt*A(1,2).*Lk);
    uk2 = M2.*(uk-3*dt*A(2,1).*Dk.*u2k1-3*dt*A(2,2).*Dk.*u2k2).*(1-dt*A(1,1).*Lk) ...
         +M2.*(uk-3*dt*A(1,1).*Dk.*u2k1-3*dt*A(1,2).*Dk.*u2k2).*(dt*A(2,1).*Lk);
    
    
    err1 = sum(abs(uk1-uk1p).^2);
    err2 = sum(abs(uk2-uk2p).^2);
    if sqrt(err1)<tol && sqrt(err2)<tol && err1+err2>=err1p+err2p
%         display(['iteration steps N = ',num2str(ll)]);
        ik=0;
        break;
    end
    
    % update states
    err1p=err1;
    err2p=err2;
    uk1p=uk1;
    uk2p=uk2;
    
    if ll==Ni
        disp(['maximum interation N = ',num2str(ll),' reached!']);
        ik=1;
    end
end