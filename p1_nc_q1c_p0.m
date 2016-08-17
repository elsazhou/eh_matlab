%P1 initialization p
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
opt7_p0=zeros(2,m);
opt7_time1=zeros(m,1);
opt7_time2=zeros(m,1);
opt7_time3=zeros(m,1);
opt7_p1=zeros(m,1);
opt7_p2=zeros(m,1);
opt8_p0=zeros(2,m);
opt8_time1=zeros(m,1);
opt8_time2=zeros(m,1);
opt8_time3=zeros(m,1);
opt8_p1=zeros(m,1);
opt8_p2=zeros(m,1);
%ini t
%load('p1_nc_q1c_ca_p025.mat');
%load('p1_nc_q1c_cb_p025.mat');
tot_iter=0;
p1s=0.005;p1e=0.3;
p2s=0.1;p2e=0.5;
for ic=1:m
    Q=Q_all(ic,:)';   
   for p01=p1s:(p1e-p1s)/25:p1e
        for p02=p2s:(p2e-p2s)/25:p2e
            p0=[p01;p02];
            [Pka1,tka1,throuput1,iter]=P1_nc_ca_p0(a,p0,Q,d,sigma,epsi);
            sumthrp1=sum(throuput1);
            if(maxthrpt7(ic)<sumthrp1)
                maxthrpt7(ic)=sumthrp1;
                opt7_p0(:,ic)=p0;
                opt7_time1(ic)=tka1(1);
                opt7_time2(ic)=tka1(2);
                opt7_time3(ic)=tka1(3);
                opt7_p1(ic)=Pka1(1);
                opt7_p2(ic)=Pka1(2);
                save p1_nc_q1c_ca_p025.mat ic Q_all opt7_time1 opt7_time2 opt7_time3 opt7_p0 maxthrpt7 opt7_p1 opt7_p2
            end
            tot_iter=tot_iter+1;
        end
    end
end

for ic=1:m
    Q=Q_all(ic,:)';
    %p1s=0.1;p1e=1;
    %p2s=0.1;p2e=1;
    for p01=p1s:(p1e-p1s)/25:p1e
        for p02=p2s:(p2e-p2s)/25:p2e
            p0=[p01;p02];
            [Pkb1,tkb1,throuput2,iter]=P1_nc_cb_p0(a,p0,Q,d,sigma,epsi);
            sumthrp2=sum(throuput2);
            if(maxthrpt8(ic)<sumthrp2)
                maxthrpt8(ic)=sumthrp2;
                opt8_p0(:,ic)=p0;
                opt8_time1(ic)=tkb1(1);
                opt8_time2(ic)=tkb1(2);
                opt8_time3(ic)=tkb1(3);
                opt8_p1(ic)=Pkb1(1);
                opt8_p2(ic)=Pkb1(2);              
                save p1_nc_q1c_cb_p025.mat ic Q_all opt8_time1 opt8_time2 opt8_time3 opt8_p0 maxthrpt8 opt8_p1 opt8_p2
            end
            tot_iter=tot_iter+1;
        end
    end
end

pli=1000*Q_all(:,1);
figure;
plot(pli,maxthrpt7,'ro-');
hold on
plot(pli,maxthrpt8,'kd-');
hold on
xlabel('Energy arrival rates of node 1(mW)')
ylabel('Troughput(bps/Hz)')
title('Comparison of Sum-Throughput (X2=100mW,ini power)')
legend('no cooperation (CA)','no cooperation (CB)')