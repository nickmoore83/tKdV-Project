% find the downstream \theta from secant method
function [theta_p,Hamil_m,Hamil_p, skew_m,skew_p,enek_m,enek_p, ik] = matching_secant_scaled(theta_m, Nw,Dm,Dp,MC)
% initialize
[~,samp_H1,samp_u,enek] = sampling_MH_scaled(Dm,Dp,Nw, theta_m, MC);
Hamil_m = mean(samp_H1);
skew_m=skewness(samp_u')';
enek_m=enek;

theta2 = -1; %theta_m;
theta1 = -2; %theta_m - 5;
[samp_H0,~,samp_u,enek] = sampling_MH_scaled(Dp,Dp,Nw, theta1, MC);
h1 = mean(samp_H0);
    
err0 = 0;
Ni=100;  % maximum interation steps
tol = 1e-1;
for ll=1:Ni
    [samp_H0,~,samp_u,enek] = sampling_MH_scaled(Dp,Dp,Nw, theta2, MC);
    h2 = mean(samp_H0);
    

    theta=theta2+(h2-Hamil_m)*(theta2-theta1)/(-h2+h1);
    err=abs(theta-theta2);
    display(['iteration i = ',num2str(ll), ', theta = ',num2str(theta), ', err = ',num2str(err)]);
    if err<tol %&& err>=err0
        display(['total iteration N = ',num2str(ll)]);
        ik=0;
        break;
    end
    h1=h2;
    theta1=theta2;
    theta2=theta;
    
    err0=err;
    if ll==Ni
        disp(['maximum interation N = ',num2str(ll),' reached!']);
        disp(['err = ',num2str(max(err)),', sum = ', num2str(sum(err>=err0))]);
        ik=1;
    end
end
theta_p = theta;
Hamil_p = h2;
skew_p=skewness(samp_u')';
enek_p=enek;