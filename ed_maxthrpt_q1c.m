%p2_1 initialization t
%Q1 is changing
clear all
clc
a=2; %path loss exponent
%Q_all=1*10^(-1).*[1 0.25;1 0.5;1 0.75;1 1;1 1.5;1 2;1 2.5;1 3;1 3.5;1 4];%x2is changing
Q_all=1*10^(-1).*[0.25 1;0.5 1;0.75 1;1 1;1.5 1;2 1;2.5 1;3 1;3.5 1;4 1];%x1 is changing
eta=0.8;
m=size(Q_all,1);
d=[1;2;1];
sigma=[10^(-4);10^(-4)];
epsi=0.0001;
%record optimal
maxthrpt1=zeros(m,1);
maxthrpt2=zeros(m,1);
opt1_t0=zeros(4,m);
opt2_t0=zeros(4,m);
opt1_time1=zeros(m,1);
opt1_time2=zeros(m,1);
opt1_time3=zeros(m,1);
opt1_time4=zeros(m,1);
opt2_time1=zeros(m,1);
opt2_time2=zeros(m,1);
opt2_time3=zeros(m,1);
opt2_time4=zeros(m,1);
opt1_rho=zeros(m,1);
opt2_rho=zeros(m,1);
%ini t
tot_iter=0;
load('ed_maxthrpt_q1c_ca.mat');
load('ed_maxthrpt_q1c_cb.mat');
for ic=m:m
    Q=Q_all(ic,:)';
    for t04=0.1:0.1:0.5
        for t01=0.2:0.1:(0.8-t04)
            for t02=0.2:0.1:(0.9-t01-t04)
                t03=1-t01-t02-t04;
                t0=[t01;t02;t03;t04];
                [process_t1,Pka1,tka1,iter1,throuput1]=ehp2_ca(a,t0,Q,eta,d,sigma,epsi);
                sumthrp1=sum(throuput1);
                if(maxthrpt1(ic)<sumthrp1)
                    maxthrpt1(ic)=sumthrp1;
                    opt1_t0(:,ic)=t0;
                    opt1_time1(ic)=tka1(1);
                    opt1_time2(ic)=tka1(2);
                    opt1_time3(ic)=tka1(3);
                    opt1_time4(ic)=tka1(4);
                    opt1_rho(ic)=Pka1(4)/Pka1(2);
                    save ed_maxthrpt_q1c_ca.mat ic Q_all opt1_time1 opt1_time2 opt1_time3 opt1_time4 opt1_t0 maxthrpt1 opt1_rho
                end              
                tot_iter=tot_iter+1;
            end
        end
    end
end
for ic=2:m
    Q=Q_all(ic,:)';
    for t04=0.05:0.1:0.2
        for t01=0.2:0.1:(0.7-t04)
            for t02=0.02:0.02:0.08
                t03=1-t01-t02-t04;
                t0=[t01;t02;t03;t04];
                [process_t2,Pkb1,tkb1,iter2,throuput2]=ehp2_cb(a,t0,Q,eta,d,sigma,epsi);
                sumthrp2=sum(throuput2);
                if(maxthrpt2(ic)<sumthrp2)
                    maxthrpt2(ic)=sumthrp2;
                    opt2_t0(:,ic)=t0;
                    opt2_time1(ic)=tkb1(1);
                    opt2_time2(ic)=tkb1(2);
                    opt2_time3(ic)=tkb1(3);
                    opt2_time4(ic)=tkb1(4);
                    opt2_rho(ic)=Pkb1(4)/Pkb1(1);
                    save ed_maxthrpt_q1c_cb.mat ic Q_all opt2_time1 opt2_time2 opt2_time3 opt2_time4 opt2_t0 maxthrpt2 opt2_rho
                end
                tot_iter=tot_iter+1;
            end
        end
    end
end

pli=1000*Q_all(:,1);
figure;
plot(pli,maxthrpt1,'ro-');
hold on
plot(pli,maxthrpt2,'md-');
hold on
xlabel('Energy arrival rates of node 1(mW)')
ylabel('Troughput(bps/Hz)')
title('Comparison of Sum-Throughput (X2=100mW,\eta=0.8)')
legend('Case A','Case B')

figure;
plot(pli,opt1_rho,'mo-');
hold on
plot(pli,opt2_rho,'bo-');
hold on
xlabel('Energy arrival rates of node 1(mW)')
ylabel('Optimal \rho')
title('X2=100mW,\eta=0.8')
legend('Case A \rhoA','Case B \rhoB')