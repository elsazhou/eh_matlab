%P1 initialization t
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
opt7_p1=zeros(m,1);
opt7_p2=zeros(m,1);
opt8_t0=zeros(3,m);
opt8_time1=zeros(m,1);
opt8_time2=zeros(m,1);
opt8_time3=zeros(m,1);
opt8_p1=zeros(m,1);
opt8_p2=zeros(m,1);
%ini t
load('p1_nc_q1c_ca_t0.mat');
%load('p1_nc_q1c_cb_t0.mat');
%step=0.005;
tot_iter=0;
for ic=9:m
    Q=Q_all(ic,:)';
    t1s=0.4;t1e=0.75;
    t2s=0.1;
   for t01=t1s:(t1e-t1s)/15:t1e
        for t02=t2s:(0.92-t01-t2s)/15:(0.9-t01)
            t03=1-t01-t02;
            t0=[t01;t02;t03];
            [Pka1,tka1,throuput1,iter]=P1_nc_ca_t0(a,t0,Q,d,sigma,epsi);
            sumthrp1=sum(throuput1);
            if(maxthrpt7(ic)<sumthrp1)
                maxthrpt7(ic)=sumthrp1;
                opt7_t0(:,ic)=t0;
                opt7_time1(ic)=tka1(1);
                opt7_time2(ic)=tka1(2);
                opt7_time3(ic)=tka1(3);
                opt7_p1(ic)=Pka1(1);
                opt7_p2(ic)=Pka1(2);
                save p1_nc_q1c_ca_t0.mat ic Q_all opt7_time1 opt7_time2 opt7_time3 opt7_t0 maxthrpt7 opt7_p1 opt7_p2
            end
            tot_iter=tot_iter+1;
        end
    end
end

for ic=1:m
    Q=Q_all(ic,:)';
    t1s=0.1;t1e=0.7;
    t2s=0.1;
   for t01=t1s:(t1e-t1s)/15:t1e
        for t02=t2s:(0.95-t01-t2s)/15:(0.95-t01)
            t03=1-t01-t02;
            t0=[t01;t02;t03];
            [Pkb1,tkb1,throuput2,iter]=P1_nc_cb_t0(a,t0,Q,d,sigma,epsi);
            sumthrp2=sum(throuput2);
            if(maxthrpt8(ic)<sumthrp2)
                maxthrpt8(ic)=sumthrp2;
                opt8_t0(:,ic)=t0;
                opt8_time1(ic)=tkb1(1);
                opt8_time2(ic)=tkb1(2);
                opt8_time3(ic)=tkb1(3);
                opt8_p1(ic)=Pkb1(1);
                opt8_p2(ic)=Pkb1(2);              
                save p1_nc_q1c_cb_t0.mat ic Q_all opt8_time1 opt8_time2 opt8_time3 opt8_t0 maxthrpt8 opt8_p1 opt8_p2
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