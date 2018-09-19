    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the source code implementation of the paper: "Robust convex
    %          approximation methods for TDOA-based localization under NLOS conditions"
    %
    % Author: Gang Wang, 
    % Ningbo University
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [y1]=sdp_ce(s,d,b,N,p1,a_)
 
    c_=[];
    C=[];
    c_conj=[];
    C_conj=[];
    k=[];
    k_conj=[];
    seita=[];
    seita_conj=[];
   v=[];
for i=1:N-1%1:N            %SDP参量 v
       v(:,i)=[zeros(1,i+1),1,zeros(1,N-i-1)].';
end
    for i=1:N-1%1:N
        c(:,i)=a_(:,i)+2*p1*v(:,i);  
        C(:,:,i)=c(:,i)*c(:,i).';
        c_conj(:,i)=a_(:,i)-2*p1*v(:,i);  
        C_conj(:,:,i)=c_conj(:,i)*c_conj(:,i).'; 
    end
  
    for i=1:N-1%1:N
         k(i)=p1^2+2*d(i)*p1+b(i);
         k_conj(i)=p1^2+2*d(i)*p1-b(i);
         seita(i)=p1^2-2*d(i)*p1+b(i);
         seita_conj(i)=p1^2-2*d(i)*p1-b(i);
    end
    for i=1:N-1%1:N
         D(:,:,i)=diag([zeros(1,i-1),-1,zeros(1,N-i-1)],0);
      
         Q1(:,:,i)=blkdiag(eye(2),D(:,:,i));
      
    end


cvx_clear
cvx_begin sdp
cvx_solver   sdpt3% sdpt3
cvx_quiet(1)
cvx_precision best %输出结果可能因此报错


variable tao(N-1,1)%tao(N-1,1)
variable y1(1+N,1)
variable Y(1+N,1+N) symmetric

minimize  (sum(tao))

subject to 


for i=1:N-1%1:N
 trace(C_conj(:,:,i)*Y)+2*k_conj(i)*c_conj(:,i).'*y1+k_conj(i)^2<=tao(i);
 
 trace(C(:,:,i)*Y)+2*seita_conj(i)*c(:,i).'*y1+seita_conj(i)^2<=tao(i);

 trace(C(:,:,i)*Y)-2*k(i)*c(:,i).'*y1+k(i)^2<=tao(i);

 trace(C_conj(:,:,i)*Y)-2*seita(i)*c_conj(:,i).'*y1+seita(i)^2<=tao(i);

end

for j=1:N-1%1:N
    trace(Q1(:,:,j)*Y)-2*s(:,j+1)'*[y1(1);y1(2)]+norm(s(:,j+1))^2==0;
    
    norm([y1(1);y1(2)]-s(:,j+1))<=y1(2+j);
   
end

[Y y1;y1.' 1]>=zeros(N+2,N+2); 

cvx_end
end
