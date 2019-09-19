theta_m =  -5:-5:-50;

Nw = [1,2,4,8,16];
MC=1e4;

for ss=1:length(Nw)
enek_m=zeros(32,length(theta_m));
enek_p=zeros(32,length(theta_m));
for ll=1:length(theta_m)
    [theta_p(ll),Hamil_m(ll),Hamil_p(ll), skew_m(:,ll),skew_p(:,ll),enek_m(:,ll),enek_p(:,ll), ik(ll)] = matching_secant_scaled(theta_m(ll), Nw(ss),1,.24,MC);
    [Nw(ss),theta_m(ll), theta_p(ll), mean(skew_m(:,ll)),mean(skew_p(:,ll))]
end

save(['matching_N',num2str(Nw(ss))]);
figure(11)
plot(theta_m,theta_p); hold on
drawnow
end
