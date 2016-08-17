%p2_1 initialization t
%Q1 is changing
%only energy cooperation
clear all
clc
a=2; %path loss exponent
%Q_all=1*10^(-1).*[1 0.25;1 0.5;1 0.75;1 1;1 1.5;1 1;1 2.5;1 3;1 3.5;1 4];%x2is changing
Q_all=1*10^(-1).*[0.25 1;0.5 1;0.75 1;1 1;1.5 1;2 1;2.5 1;3 1;3.5 1;4 1];%x1 is changing
eta=0.8;
m=size(Q_all,1);
d=[1;2;1];
sigma=[10^(-4);10^(-4)];
epsi=0.0001;
%record optimal
maxthrpt5=zeros(m,1);
maxthrpt6=zeros(m,1);
opt5_t0=zeros(3,m);
opt5_time1=zeros(m,1);
opt5_time2=zeros(m,1);
opt5_time3=zeros(m,1);
opt6_t0=zeros(3,m);
opt6_time1=zeros(m,1);
opt6_time2=zeros(m,1);
opt6_time3=zeros(m,1);
%ini t

tot_iter=0;
load('e_maxthrpt_q1c_ca.mat');
load('e_maxthrpt_q1c_cb.mat');
step=0.005;
for ic=m:m
    Q=Q_all(ic,:)';
    for t01=opt5_time1(ic)-opt5_time3(ic)/2:step:opt5_time1(ic)+opt5_time3(ic)/2-step
        for t02=opt5_time2(ic)-opt5_time3(ic)/2:step:opt5_time2(ic)+opt5_time3(ic)/2-step
            t03=1-t01-t02;
            t0=[t01;t02;t03];
            [Pka1,tka1,throuput1]=ehp2_ec_ca(a,t0,Q,d,sigma,epsi,eta);
            sumthrp1=sum(throuput1);
            if(maxthrpt5(ic)<sumthrp1)
                maxthrpt5(ic)=sumthrp1;
                opt5_t0(:,ic)=t0;
                opt5_time1(ic)=tka1(1);
                opt5_time2(ic)=tka1(2);
                opt5_time3(ic)=tka1(3);
                save e_maxthrpt_q1c_ca.mat ic Q_all opt5_time1 opt5_time2 opt5_time3 opt5_t0 maxthrpt5 step
            end
            tot_iter=tot_iter+1;
        end
    end
end
for ic=7:m
    Q=Q_all(ic,:)';
    for t01=(opt6_time1(ic)-opt6_time3(ic)/2):step:(opt6_time1(ic)+opt6_time3(ic)/2-0.01)
        for t03=(opt6_time2(ic)-opt6_time3(ic)/2):step:(opt6_time3(ic)+opt6_time3(ic)/2-0.01)
            t02=1-t01-t03;
            t0=[t01;t02;t03];
            [Pkb1,tkb1,throuput2]=ehp2_ec_cb(a,t0,Q,d,sigma,epsi,eta);
            sumthrp2=sum(throuput2);
            if(maxthrpt6(ic)<sumthrp2)
                maxthrpt6(ic)=sumthrp2;
                opt6_t0(:,ic)=t0;
                opt6_time1(ic)=tkb1(1);
                opt6_time2(ic)=tkb1(2);
                opt6_time3(ic)=tkb1(3);
                save e_maxthrpt_q1c_cb.mat ic Q_all opt6_time1 opt6_time2 opt6_time3 opt6_t0 maxthrpt6 step
            end
            tot_iter=tot_iter+1;
        end
    end
end
pli=1000*Q_all(:,1);
figure;
plot(pli,maxthrpt5,'ro-');
hold on
plot(pli,maxthrpt6,'md-');
hold on
xlabel('Energy arrival rates of node 1(mW)')
ylabel('Troughput(bps/Hz)')
title('Comparison of Sum-Throughput (X2=100mW)')
legend('only energy cooperation (CA)','only energy cooperation (CB)')