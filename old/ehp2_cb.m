%%%%%Problem 2 Max throughput
% parameters setting
%     h -- channel gains,h=10^-3*d^-a;  h1<h2<h3
% sigma -- noise variance, 
%   P_c -- circuit power
% eta   -- energy tranfer efficiency
% Q=[Q1,Q2] -- average energy havesting rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [process_t,power,time,iter,throuput]=ehp2_cb(a,t0,Q,eta,d,sigma,epsi)
h=d.^(-a);
% Alternating Convex Programming
% variables t(4) P(3)
tk1=t0;
tk2=t0;
T=sum(t0);
tk=zeros(4,1);
Pk2=ones(4,1);
Pk=zeros(4,1);
iter=0;
start_t=cputime;
while((norm(Pk-Pk2) > epsi || norm(tk-tk2) > epsi)&& iter < 15)
    tk=tk2;
    Pk=Pk2;
    [thr1,Pk2]=ehp2_fixtB(tk2,h,sigma,eta,Q);
    [thr2,tk2]=ehp2_fixpB(T,Pk2,h,sigma,eta,Q); %cvx solutions
    iter=iter+1;
end
end_t=cputime;
process_t=end_t-start_t;
throuput=thr2;
power=Pk2;
time=tk2;