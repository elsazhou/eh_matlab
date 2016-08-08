% Optimization Project
% Problem 2
% if time is fixed, find an optimal Power
% input variables
%     h -- channel gains
% sigma -- noise variance
% eta   -- energy tranfer efficiency
% Q1 -- average energy havesting rate of node 1
function [thrpt,power]=ehp2_fixtndcA(t,h,sigma,Q)
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
    variable P(3);
    variable B;
    minimize(-B-t1*log(1+P(1)*r1));
    subject to 
        t2*log(1+P(2)*r2)+t3*log(1+P(3)*r1)>=B;
        t2*log(1+P(2)*r3)>=B;
        Q1*t4>=P(1)*t1;
        Q2*(t1+t4)>=P(2)*t2;
        Q1*(t1+t2+t4)-P(1)*t1>=P(3)*t3;
        P(1)>= 0;
        P(2)>= 0;     
        P(3)>= 0;
cvx_end
power=P;
thrpt=zeros(2,1);
thrpt(1)=t1*log(1+P(1)*r1);
thrpt(2)=B;