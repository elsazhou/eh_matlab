function [Pk1,tk,throuput]=ehp2_nc_ca(a,t0,Q,d,sigma,epsi)
r1=d(1)^(-a)/sigma(1);
r2=d(2)^(-a)/sigma(1);
r=[r1;r2];
tk1=t0;
tk=zeros(3,1);
Pk1=ones(2,1);
Pk=zeros(2,1);
T=1;
iter=0;
while((norm(Pk-Pk1) > epsi || norm(tk-tk1) > epsi)&&iter<20)
    tk=tk1;
    Pk=Pk1;
    [thr1,Pk1]=ehp2_nc_fixtA(tk1,r,Q);
    [thr2,tk1]=ehp2_nc_fixpA(T,Pk1,r,Q); %cvx solutions
    iter=iter+1;
end
throuput=thr2;
tk=tk1;