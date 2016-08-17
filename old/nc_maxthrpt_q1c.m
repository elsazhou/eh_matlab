%p2_1 initialization t
%Q2 is changing
%no data relay no energy
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
maxthrpt7=zeros(m,1);
maxthrpt8=zeros(m,1);
opt7_t0=zeros(3,m);
opt7_time1=zeros(m,1);
opt7_time2=zeros(m,1);
opt7_time3=zeros(m,1);
opt8_t0=zeros(3,m);
opt8_time1=zeros(m,1);
opt8_time2=zeros(m,1);
opt8_time3=zeros(m,1);
%ini t
load('nc_maxthrpt_q1c_ca.mat');
load('nc_maxthrpt_q1c_cb.mat');
step=0.005;
tot_iter=0;
for ic=9:m
    Q=Q_all(ic,:)';
   for t01=(opt7_time1(ic)-0.1):step:(opt7_time1(ic)+opt7_time3(ic)/2-0.01)
        for t03=0.005:step:(opt7_time3(ic)+0.1)
            t02=1-t01-t03;
            t0=[t01;t02;t03];
            [Pka1,tka1,throuput1]=ehp2_nc_ca(a,t0,Q,d,sigma,epsi);
            sumthrp1=sum(throuput1);
            if(maxthrpt7(ic)<sumthrp1)
                maxthrpt7(ic)=sumthrp1;
                opt7_t0(:,ic)=t0;
                opt7_time1(ic)=tka1(1);
                opt7_time2(ic)=tka1(2);
                opt7_time3(ic)=tka1(3);
                save nc_maxthrpt_q1c_ca.mat ic Q_all opt7_time1 opt7_time2 opt7_time3 opt7_t0 maxthrpt7 step
            end
            tot_iter=tot_iter+1;
        end
    end
end

for ic=1:m
    Q=Q_all(ic,:)';
    for t01=(opt8_time1(ic)-0.1):step:(opt8_time1(ic)+opt8_time3(ic)/2-0.01)
        for t03=0.005:step:(opt8_time3(ic)+0.1)
            t02=1-t01-t03;
            t0=[t01;t02;t03];
            [Pkb1,tkb1,throuput2]=ehp2_nc_cb(a,t0,Q,d,sigma,epsi);
            sumthrp2=sum(throuput2);
            if(maxthrpt8(ic)<sumthrp2)
                maxthrpt8(ic)=sumthrp2;
                opt8_t0(:,ic)=t0;
                opt8_time1(ic)=tkb1(1);
                opt8_time2(ic)=tkb1(2);
                opt8_time3(ic)=tkb1(3);
                save nc_maxthrpt_q1c_cb.mat ic Q_all opt8_time1 opt8_time2 opt8_time3 opt8_t0 maxthrpt8 step
            end
            tot_iter=tot_iter+1;
        end
    end
end

pli=1000*Q_all(:,1);
figure;
plot(pli,maxthrpt7,'ro-');
hold on
plot(pli,maxthrpt8,'md-');
hold on
xlabel('Energy arrival rates of node 1(mW)')
ylabel('Troughput(bps/Hz)')
title('Comparison of Sum-Throughput (X2=100mW)')
legend('no cooperation (CA)','no cooperation (CB)')