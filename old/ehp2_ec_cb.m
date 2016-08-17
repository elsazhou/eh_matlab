function [Pk1,tk1,throuput]=ehp2_ec_cb(a,t0,Q,d,sigma,epsi,eta)
h=d.^(-a);
% Alternating Convex Programming
tk1=t0;
tk=zeros(3,1);
Pk1=zeros(2,1);
Pk=ones(2,1);
iter=0;
T=sum(t0);
while((norm(Pk-Pk1) > epsi || norm(tk-tk1) > epsi)&& iter<15)
    tk=tk1;
    Pk=Pk1;
    [thr1,Pk1]=ehp2_fixt_ecB(tk1,h,sigma,Q,eta);
    [thr2,tk1]=ehp2_fixp_ecB(T,Pk1,h,sigma,Q,eta);
    iter=iter+1;
end
throuput=thr2;