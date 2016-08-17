function [Pk1,tk,throuput,iter]=P1_nc_ca_t0(a,t0,Q,d,sigma,epsi)
r1=d(1)^(-a)/sigma(1);
r2=d(2)^(-a)/sigma(1);
r=[r1;r2];
tk1=t0;
tk=zeros(3,1);
Pk1=ones(2,1);
Pk=zeros(2,1);
T=1;
iter=0;
while((norm(Pk-Pk1) > epsi || norm(tk-tk1) > epsi))
    tk=tk1;
    Pk=Pk1;
    [thr1,Pk1]=P1_nc_fixtA(tk1,r,Q);
    [thr2,tk1]=P1_nc_fixpA(Pk1,r,Q); %sedumi
    iter=iter+1;
end
throuput=thr2;
tk=tk1;