% Optimization Project
% Problem 2
% if power is fixed, find an optimal time using CVX
% input variables
%     h -- channel gains
% sigma -- noise variance
% eta   -- energy tranfer efficiency
% Q1 -- average energy havesting rate of node 1
function [thrp,power]=ehp2_fixtB(t,h,sigma,eta,Q)
t1=t(1);
t2=t(2);
t3=t(3);
t4=t(4);
r1=h(1)/sigma(1);
r2=h(2)/sigma(1);
r3=h(3)/sigma(2);
Q1=Q(1);
Q2=Q(2);
cvx_begin quiet
    variable P(4);
    variable B;
    minimize(-B-t3*log(1+P(3)*r1));
    subject to
        t1*log(1+P(1)*r2)+t2*log(1+P(2)*r1)>=B;
        t1*log(1+(P(1)-P(4))*r3)>=B;
        Q2*t4>=P(1)*t1;
        Q1*t4+(eta*h(3)*P(4)+Q1)*t1>=P(2)*t2;
        Q1*t4-P(2)*t2+Q1*(t1+t2)+eta*h(3)*P(4)*t1>=P(3)*t3;      
        P(1)>=0;
        P(2)>=0;
        P(3)>=0;
        P(4)>=0;
        P(1)>=P(4);
cvx_end
power=P;
thrp=zeros(2,1);
thrp(2)=B;
thrp(1)=t3*log(1+P(3)*r1);