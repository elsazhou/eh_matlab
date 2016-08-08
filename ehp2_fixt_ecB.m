function [thrp,power]=ehp2_fixt_ecB(t,h,sigma,Q,eta)
t1=t(1);
t2=t(2);
t3=t(3);
r1=h(1)/sigma(1);
r2=h(2)/sigma(1);
Q1=Q(1);
Q2=Q(2);
cvx_begin quiet
    variable P(2);
    minimize(-t1*log(1+P(1)*r2)-t2*log(1+P(2)*r1));
    subject to 
        Q2*t3>=P(1)*t1;
        Q1*t3+(eta*h(3)*P(1)+Q1)*t1>=P(2)*t2;
        P(1)>=0;
        P(2)>=0;
cvx_end
power=P;
thrp=zeros(2,1);
thrp(1)=t2*log(1+P(2)*r1);
thrp(2)=t1*log(1+P(1)*r2);