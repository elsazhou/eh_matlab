%p2_1 initialization t
%Q2 is changing
%only data relay
clear all
clc
a=2; %path loss exponent
Q_all=1*10^(-1).*[1 0.25;1 0.5;1 0.75;1 1;1 1.5;1 1;1 2.5;1 3;1 3.5;1 4];%x2is changing
%Q_all=1*10^(-1).*[0.25 1;0.5 1;0.75 1;1 1;1.5 1;2 1;2.5 1;3 1;3.5 1;4 1];%x1 is changing
eta=0.8;
m=size(Q_all,1);
d=[1;2;1];
sigma=[10^(-4);10^(-4)];
epsi=0.0001;
%record optimal
maxthrpt3=zeros(m,1);
maxthrpt4=zeros(m,1);
opt3_t0=zeros(4,m);
opt3_time1=zeros(m,1);
opt3_time2=zeros(m,1);
opt3_time3=zeros(m,1);
opt3_time4=zeros(m,1);
opt4_t0=zeros(4,m);
opt4_time1=zeros(m,1);
opt4_time2=zeros(m,1);
opt4_time3=zeros(m,1);
opt4_time4=zeros(m,1);
%ini t
tot_iter=0;
step=0.005;
load('d_maxthrpt_q2c_ca.mat');
load('d_maxthrpt_q2c_cb.mat');

for ic=1:m
    Q=Q_all(ic,:)';
    for t04=(opt3_time4(ic)-0.1):step:(opt3_time4(ic)+0.1)
        for t02=(opt3_time2(ic)-0.1):step:(opt3_time2(ic)+0.1)
            for t03=0.03:step:0.02
                t01=1-t03-t02-t04;
                t0=[t01;t02;t03;t04];
                [process_t1,Pka1,tka1,iter1,throuput1]=ehp2_ndc_ca(a,t0,Q,d,sigma,epsi);
                sumthrp1=sum(throuput1);
                if(maxthrpt3(ic)<sumthrp1)
                    maxthrpt3(ic)=sumthrp1;
                    opt3_t0(:,ic)=t0;
                    opt3_time1(ic)=tka1(1);
                    opt3_time2(ic)=tka1(2);
                    opt3_time3(ic)=tka1(3);
                    opt3_time4(ic)=tka1(4);
                    save d_maxthrpt_q2c_ca.mat ic Q_all opt3_time1 opt3_time2 opt3_time3 opt3_time4 opt3_t0 maxthrpt3 step
                end
                tot_iter=tot_iter+1;
            end
        end
    end
end
for ic=1:m
    Q=Q_all(ic,:)';
    for t04=(opt4_time4(ic)-0.1):step:(opt4_time4(ic)+0.1)
        for t01=(opt4_time1(ic)-0.1):step:(opt4_time1(ic)+0.1)
            for t02=0.03:step:0.02              
                t03=1-t01-t02-t04;
                t0=[t01;t02;t03;t04];                
                [process_t2,Pkb1,tkb1,iter2,throuput2]=ehp2_ndc_cb(a,t0,Q,d,sigma,epsi);
                sumthrp2=sum(throuput2);
                if(maxthrpt4(ic)<sumthrp2)
                    maxthrpt4(ic)=sumthrp2;
                    opt4_t0(:,ic)=t0;
                    opt4_time1(ic)=tkb1(1);
                    opt4_time2(ic)=tkb1(2);
                    opt4_time3(ic)=tkb1(3);
                    opt4_time4(ic)=tkb1(4);
                    save d_maxthrpt_q2c_cb.mat ic Q_all opt4_time1 opt4_time2 opt4_time3 opt4_time4 opt4_t0 maxthrpt4 step
                end
                tot_iter=tot_iter+1;
            end
        end
    end
end

pli=1000*Q_all(:,2);
figure;
plot(pli,maxthrpt3,'ro-');
hold on
plot(pli,maxthrpt4,'md-');
hold on
xlabel('Energy arrival rates of node 2(mW)')
ylabel('Troughput(bps/Hz)')
title('Comparison of Sum-Throughput (X1=100mW)')
legend('only data relay (CA)','only data relay (CB)')