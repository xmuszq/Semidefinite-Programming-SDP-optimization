clc
clear all
close all
tic


MC=10; % how many runs
p=10; % maximal NLOS value
N=8; % number of anchors
t0=0.1;%���Ȳ�������t1һ��ֵ��֤t0����Ϊ����

c=1;

%sigma=0.3;
%x_(1)=0;
%y_(1)=0;
%s=[];
%s(:,1)=[x_(1);y_(1)];%�ο�ê�ڵ�����

s(:,1)=[-20;-20];%anchors
s(:,2)=[-20;20];
s(:,3)=[20;-20];
s(:,4)=[20;20];
s(:,5)=[-20;0];
s(:,6)=[0;-20];
s(:,7)=[0;20];
s(:,8)=[20;0];


cond=[1,5,6,8] %number of NLOS
% loop for nlos for 1~5 and 6~8 respectively
for ci=1:2
    if ci==1
        nrange=[1,5]
    else
        nrange=[6,8]
    end
    
    
    K=0.1:0.1:1;%�������������
    for m=1:10
        sigma=0.1*m;
        %Q=0.5*sigma^2*(eye(N-1)+ones(N-1,N-1));%����Э�������
        Q=0.5*sigma^2*(eye(N)+ones(N,N));%����Э�������
        parfor j=1:MC
            j;
            
            NLOS_number = randi([6,8],1,1);  %randsrc(1,1,randperm(5));%��ѡ·��
            %n=gauss_samples(zeros(N-1,1),Q,1); % get the noise
            n=gauss_samples(zeros(N,1),Q,1); % get the noise
            
            x=unifrnd(-25,25);%Ŀ��ڵ�
            y=unifrnd(-25,25);
            X=[x;y];
            
            w = [];
            for i=1:N%�������Ӿ���� get the nlos bias for 8 links
                w(i)=unifrnd(0,p);
            end
            
            g_index=randperm(8,NLOS_number); % assign nlos bias to certain links
            g = zeros(N,1);
            g(g_index) = 1;
            
            
            %%%%%%%%%%%% SDP wang 2016 %%%%%%%%%%%%%%
            e=[];
            d=[];
            b=[];
            for i=1:N-1%1:N;%NLOS������ֵ
                e(i,1)=g(i+1)*w(i+1)-g(1)*w(1);%w(i+1)-w(1)
                % e_(i-1,1)=w_(i)-w_(1);
                d(i)=norm(X-s(:,i+1))-norm(X-s(:,1))+n(i)+e(i);
                %d(i)=norm(X-s(:,i))-norm(X-s0)+n(i)+e(i);
                b(i)=-d(i)^2-norm(s(:,i+1))^2+norm(s(:,1))^2;
                %b(i-1)=-d(i-1)^2-norm(s(:,i))^2+norm(s0)^2;
            end
            
            a_=[];
            for i=1:N-1%2:N
                %a_(:,i-1)=[2*(s(:,1)-s(:,i))',zeros(1,i-1),-2*d(i-1),zeros(1,N-i)].';
                a_(:,i)=[2*(s(:,1)-s(:,i+1))',zeros(1,i-1),-2*d(i),zeros(1,N-i-1)].';
            end
            y0=sdp_ce(s,d,b,N,p,a_);
            X1=[y0(1);y0(2)];
            t00(j, 1)=norm(X1-X);
            
            
            
            %%%%%%%%%%%% SDP CL chen %%%%%%%%%%%%%%%%%
            e = [];
            d = [];
            t1= [];
            for i=2:N%NLOS��TDOA���� generate TDOA with NLosֵ
                e(i-1,1)=g(i)*w(i)-g(1)*w(1);%w(i)-w(1);
                % e_(i-1,1)=w_(i)-w_(1);
                d(i-1)=norm(X-s(:,i))-norm(X-s(:,1))+n(i-1)+e(i-1);
                % t1(i)=d(i-1)/(3*10^8)+t1(1);
            end
            %t2=t1*10^8;
            %t1(1)=norm(X-s(:,1))/c+t0;
            t1(1)=norm(X-s(:,1))/c+t0;% ???????it's better not to add e(1) n(1)
            %t1(1)=abs(min(d))+15;
            %t1(1)=0.1;
            
            t1(2:N)=d/c+t1(1);
            t2=t1;
            %  t2
            y1=SDP_CL_chen(c,s,N,t2);
            t11(j,1)=norm(y1-X);
            
            
            
            %%%%%%%%%%%%%%% SDP CL SU %%%%%%%%%%%%%%%%%%
            % generate TOA
            toa=zeros(N, 1);
            for i=1:N
                toa(i)= norm(X-s(:,i))+g(i)*w(i)+n(i);
            end
            Cp=2.99792458 ;
            y2=SDP_CL_SU(toa, s, t0, Cp);
            t22(j,1)=norm(y2-X);
            
            
            
            %%%%%%%%%%%%%%% SDP  %%%%%%%%%%%%%%%%%%
            c_speed=1%2.99792458; %m/s
            y3=SDP(toa, s, c_speed)
            t33(j,1)=norm(y3-X);                 
            
            
            
        end
        T0(m)=sum(t00)/MC;
        T1(m)=sum(t11)/MC;
        T2(m)=sum(t22)/MC;
        T3(m)=sum(t33)/MC;
        
    end
    
    
    
    
    if ci==1
        save res_me-chen_nlos1to5 T1
        save res_me-su_nlos1to5 T2
        save res_wang_nlos1to5 T0
        save res_liu_nlos1to5 T3

        figure(1)
        plot(K, T0, 'b-')
        hold on;
        plot(K, T1,'g-')
        plot(K, T2, 'r-')
        plot(K, T3, 'c-')
        title('nlos1 1 - 5')
        xlabel('sigma');
        ylabel('RMSE');
        legend('SDP-robust wang', 'Proposed-code from Chen', 'Proposed', 'SDP - Wang')
        
    else
        save res_me-chen_nlos6to8 T1
        save res_me-su_nlos6to8 T2
        save res_wang_nlos6to8 T0
        save res_liu_nlos6to8 T3
        
        figure(2)
        plot(K, T0, 'b-')
        hold on;
        plot(K,T1,'g-')
        plot(K, T2, 'r-')
        plot(K, T3, 'c-')
        
        title('nlos1 6 -8')
        xlabel('sigma');
        ylabel('RMSE');
        legend('SDP-robust wang', 'Proposed-code from Chen', 'Proposed - code from Su', 'SDP - Wang')
        
    end
grid on
end
toc