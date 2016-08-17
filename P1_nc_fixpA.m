function [thr2,tk1]=P1_nc_fixpA(p,r,Q) %sedumi
P1=p(1);
P2=p(2);
r1=r(1);
r2=r(2);
Q1=Q(1);
Q2=Q(2);
thr2=zeros(2,1);tk1=zeros(3,1);
b=[log(1+P1*r1);log(1+P2*r2)];
A=[P1+Q1,Q1;0,P2+Q2;-1 0;0 -1;1 1];
c=[Q1;Q2;0;0;1];
[x,y]=sedumi(A,b,c);
tk1(1)=y(1);
tk1(2)=y(2);
tk1(3)=1-tk1(1)-tk1(2);
thr2(1)=tk1(1)*log(1+P1*r1);
thr2(2)=tk1(2)*log(1+P2*r2);