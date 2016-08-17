% parameters setting
%     h -- channel gains,h=10^-3*d^-a;  h1<h2<h3
% sigma -- noise variance, 
% eta   -- energy tranfer efficiency
% Q1 -- average energy havesting rate of node 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [process_t,P,t,iter,throuput]=ehp2_ndc_cb(a,t0,Q,d,sigma,epsi)
 %distance
h=10^(0).*d.^(-a);
% Alternating Convex Programming
% variables t(4) P(3)
%tk1=ones(3,1);
tk1=t0;
tk=zeros(4,1);
Pk1=zeros(3,1);
Pk=ones(3,1);
iter=0;
T=sum(t0);
start_time=cputime;
while((norm(Pk-Pk1) > epsi || norm(tk-tk1) > epsi) && iter<15)
    tk=tk1;
    Pk=Pk1;
    [thr1,Pk1]=ehp2_fixtndcB(tk1,h,sigma,Q);
    [thr2,tk1]=ehp2_fixpndcB(T,Pk1,h,sigma,Q);
    iter=iter+1;
end
end_time=cputime;
process_t=end_time-start_time;
P=Pk1;
t=tk1;
throuput=thr2;