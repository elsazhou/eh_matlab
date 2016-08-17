% non cooperation
% Problem 2
% if time is fixed, find an optimal Power using cvx
% input variables
%     h -- channel gains
% sigma -- noise variance
% Q -- average energy havesting rates of nodes
function [thrpt,power]=ehp2_nc_fixtA(t,r,Q)
t1=t(1);
t2=t(2);
t3=t(3);
r1=r(1);
r2=r(2);
Q1=Q(1);
Q2=Q(2);
cvx_begin quiet
    variable P(2);
    minimize(-t1*log(1+P(1)*r1)-t2*log(1+P(2)*r2));
    subject to 
        Q1*t3>=P(1)*t1;
        Q2*(t3+t1)>=P(2)*t2;          
        P(1)>=0;
        P(2)>=0;
cvx_end
power=P;
thrpt=zeros(2,1);
thrpt(1)=t1*log(1+P(1)*r1);
thrpt(2)=t2*log(1+P(2)*r2);