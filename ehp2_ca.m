%%%%%Problem 2 Max throughput
% parameters setting
%     h -- channel gains,h=d^-a;  
% sigma -- noise variance, 
%   P_c -- circuit power
% eta   -- energy tranfer efficiency
% Q=[Q1,Q2] -- average energy havesting rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [process_t,Pk1,tk1,iter,throuput]=ehp2_ca(a,t0,Q,eta,d,sigma,epsi)
h=d.^(-a);
% Alternating Convex Programming
% variables t(4) P(3)
tk1=t0;
tk2=t0;
T=sum(t0);
tk=zeros(4,1);
Pk1=ones(4,1);
Pk=zeros(4,1);
iter=0;
start_t=cputime;
while((norm(Pk-Pk1) > epsi || norm(tk-tk1) > epsi) && iter < 25)
    tk=tk1;
    Pk=Pk1;
    [thr1,Pk1]=ehp2_fixtA(tk1,h,sigma,eta,Q);
    [thr2,tk1]=ehp2_fixpA(T,Pk1,h,sigma,eta,Q);
    iter=iter+1;
end
end_t=cputime;
process_t=end_t-start_t;
throuput=thr2;