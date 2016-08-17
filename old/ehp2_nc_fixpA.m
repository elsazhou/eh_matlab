% non cooperation
% Problem 2
% if time is fixed, find an optimal Power using cvx
% input variables
%     h -- channel gains
% sigma -- noise variance
% Q -- average energy havesting rates of nodes
function [thrp,time]=ehp2_nc_fixpA(T,p,r,Q)
P1=p(1);
P2=p(2);
r1=r(1);
r2=r(2);
Q1=Q(1);
Q2=Q(2);
thrp=zeros(2,1);
cvx_begin quiet
    variable t(3);
    minimize(-t(1)*log(1+P1*r1)-t(2)*log(1+P2*r2));
    subject to
        Q1*t(3)>=P1*t(1);
        Q2*(t(3)+t(1))>=P2*t(2);          
        t(1)+t(2)+t(3)==T;
        t(1)>= 0;
        t(2)>= 0;
        t(3)>= 0;
cvx_end
time=t;
thrp(1)=t(1)*log(1+P1*r1);
thrp(2)=t(2)*log(1+P2*r2);