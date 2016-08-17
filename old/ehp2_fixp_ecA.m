function [thrp,time]=ehp2_fixp_ecA(T,p,h,sigma,Q,eta)
P1=p(1);
P2=p(2);
r1=h(1)/sigma(1);
r2=h(2)/sigma(1);
Q1=Q(1);
Q2=Q(2);
cvx_begin quiet
    variable t(3);
    minimize(-t(1)*log(1+P1*r1)-t(2)*log(1+P2*r2));
    subject to 
        Q1*t(3)>=P1*t(1);
        Q2*t(3)+(eta*h(3)*P1+Q2)*t(1)>=P2*t(2);
        t(1)>=0;
        t(2)>=0;
        t(3)>=0;
        t(1)+t(2)+t(3)==T;
cvx_end
time=t;
thrp=zeros(2,1);
thrp(2)=t(2)*log(1+P2*r2);
thrp(1)=t(1)*log(1+P1*r1);