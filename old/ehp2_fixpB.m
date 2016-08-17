% Optimization Project
% Problem 2
% if power is fixed, find an optimal time using CVX
% input variables
%     h -- channel gains
% sigma -- noise variance
% eta   -- energy tranfer efficiency
% Q1 -- average energy havesting rate of node 1
function [thrp,time]=ehp2_fixpB(T,p,h,sigma,eta,Q)
P1=p(1);
P2=p(2);
P3=p(3);
P4=p(4);
r1=h(1)/sigma(1);
r2=h(2)/sigma(1);
r3=h(3)/sigma(2);
Q1=Q(1);
Q2=Q(2);
thrp=zeros(2,1);
cvx_begin quiet
    variable t(4);
    variable B;
    minimize(-B-t(3)*log(1+P3*r1));
    subject to
        t(1)*log(1+P1*r2)+t(2)*log(1+P2*r1)>=B;
        t(1)*log(1+(P1-P4)*r3)>=B;
        Q2*t(4)>=P1*t(1);
        Q1*t(4)+(eta*h(3)*P4+Q1)*t(1)>=P2*t(2);
        Q1*t(4)-P2*t(2)+Q1*(t(1)+t(2))+eta*h(3)*P4*t(1)>=P3*t(3);      
        t(1)+t(2)+t(3)+t(4) == T;
        t(1)>= 0;
        t(2)>= 0;
        t(3)>= 0;
        t(4)>= 0;
cvx_end
time=t;
thrp(2)=B;
thrp(1)=t(3)*log(1+P3*r1);