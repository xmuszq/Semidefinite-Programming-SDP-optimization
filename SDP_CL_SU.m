function [x]=SDP_CL_SU_1(d, anchor_posi, dt, Cp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % d - the real TOA including the error and NLOS bias
    % dt - a small positive value to compensate the timing
    % anchor_posi - the anchor positions (hard-coded for 8 anchors)
    % x - the estimated position
    % Cp - the speed of light and is default to 2.99792
    
    % Implementation of the paper:"Semidefinite Programming for NLOS Error
    %       Mitigation in TDOA Localization"
    % Author: Zhenqiang Su
    % Email: xiamenhai@gmail.com
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %speed of light
    Cp=Cp;%2.99792;
    %general parameter for algorithms
    v=[1 1 1 1 1 1 1 1];
    u=[.2 .2 .2 .2 .2 .2 .2 .2];% u>2*variance
    u=u*1.5;
    umuso=0.00; umuso_z=0.00; extra=0.0;

    d1_observe = d(1);
    d2_observe = d(2);
    d3_observe = d(3);
    d4_observe = d(4);
    d5_observe = d(5);
    d6_observe = d(6);
    d7_observe = d(7);
    d8_observe = d(8);

    % TOA distance to timing
    t1=	d1_observe/Cp;
    t2=	d2_observe/Cp;
    t3=	d3_observe/Cp;
    t4=	d4_observe/Cp;
    t5=	d5_observe/Cp;
    t6=	d6_observe/Cp;
    t7=	d7_observe/Cp;
    t8=	d8_observe/Cp;
    tmin=min([t1,t2,t3,t4,t5,t6,t7,t8]);

    % remove the TOA information and add a constant positive dt
    % after this operation only TDOA information is reserved
    t1=	t1-tmin+dt;
    t2=	t2-tmin+dt;
    t3=	t3-tmin+dt;
    t4=	t4-tmin+dt;
    t5=	t5-tmin+dt;
    t6=	t6-tmin+dt;
    t7=	t7-tmin+dt;
    t8=	t8-tmin+dt;
    
    index=find([t1,t2,t3,t4,t5,t6,t7,t8]==tmin);
    alf_sdp_rev=5; alf=0.01;%0.01 is better than 0.1 and better than 0.02
    %%%%%%%%%%%%%%%%%%%%%%%% CpVX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin sdp quiet
    %cvx_solver sedumi
    variable x(2,1)
    variable c(8,1)
    variable z
    variable h(8,1)
    variable t0
    variable t0_s
    variable f

    minimize  v(1)*(Cp^2*t0_s-2*Cp^2*t1*t0+Cp^2*t1^2-h(1)-c(1))^2+alf*c(1)^2...
            +v(2)*(Cp^2*t0_s-2*Cp^2*t2*t0+Cp^2*t2^2-h(2)-c(2))^2+alf*c(2)^2...
            +v(3)*(Cp^2*t0_s-2*Cp^2*t3*t0+Cp^2*t3^2-h(3)-c(3))^2+alf*c(3)^2...
            +v(4)*(Cp^2*t0_s-2*Cp^2*t4*t0+Cp^2*t4^2-h(4)-c(4))^2+alf*c(4)^2...
            +v(5)*(Cp^2*t0_s-2*Cp^2*t5*t0+Cp^2*t5^2-h(5)-c(5))^2+alf*c(5)^2...
            +v(6)*(Cp^2*t0_s-2*Cp^2*t6*t0+Cp^2*t6^2-h(6)-c(6))^2+alf*c(6)^2...
            +v(7)*(Cp^2*t0_s-2*Cp^2*t7*t0+Cp^2*t7^2-h(7)-c(7))^2+alf*c(7)^2...
            +v(8)*(Cp^2*t0_s-2*Cp^2*t8*t0+Cp^2*t8^2-h(8)-c(8))^2+alf*c(8)^2...
            +alf*(t0_s*8*8)^2 +f^2;

    subject to

    h(1)==[anchor_posi(1,1)  anchor_posi(2,1)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,1); anchor_posi(2,1);-1];
    h(2)==[anchor_posi(1,2)  anchor_posi(2,2)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,2); anchor_posi(2,2);-1];
    h(3)==[anchor_posi(1,3)  anchor_posi(2,3)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,3); anchor_posi(2,3);-1];
    h(4)==[anchor_posi(1,4)  anchor_posi(2,4)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,4); anchor_posi(2,4);-1];
    h(5)==[anchor_posi(1,5)  anchor_posi(2,5)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,5); anchor_posi(2,5);-1];
    h(6)==[anchor_posi(1,6)  anchor_posi(2,6)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,6); anchor_posi(2,6);-1];
    h(7)==[anchor_posi(1,7)  anchor_posi(2,7)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,7); anchor_posi(2,7);-1];
    h(8)==[anchor_posi(1,8)  anchor_posi(2,8)  -1]*[ 1 0 x(1);0 1 x(2);x(1) x(2) z]*[anchor_posi(1,8); anchor_posi(2,8);-1];

    [1 0 anchor_posi(1,1)-x(1);0 1 anchor_posi(2,1)-x(2);anchor_posi(1,1)-x(1) anchor_posi(2,1)-x(2) Cp^2*t0_s-2*Cp^2*t1*t0+Cp^2*t1^2+f]>=0;
    [1 0 anchor_posi(1,2)-x(1);0 1 anchor_posi(2,2)-x(2);anchor_posi(1,2)-x(1) anchor_posi(2,2)-x(2) Cp^2*t0_s-2*Cp^2*t2*t0+Cp^2*t2^2+f]>=0;
    [1 0 anchor_posi(1,3)-x(1);0 1 anchor_posi(2,3)-x(2);anchor_posi(1,3)-x(1) anchor_posi(2,3)-x(2) Cp^2*t0_s-2*Cp^2*t3*t0+Cp^2*t3^2+f]>=0;
    [1 0 anchor_posi(1,4)-x(1);0 1 anchor_posi(2,4)-x(2);anchor_posi(1,4)-x(1) anchor_posi(2,4)-x(2) Cp^2*t0_s-2*Cp^2*t4*t0+Cp^2*t4^2+f]>=0;
    [1 0 anchor_posi(1,5)-x(1);0 1 anchor_posi(2,5)-x(2);anchor_posi(1,5)-x(1) anchor_posi(2,5)-x(2) Cp^2*t0_s-2*Cp^2*t5*t0+Cp^2*t5^2+f]>=0;
    [1 0 anchor_posi(1,6)-x(1);0 1 anchor_posi(2,6)-x(2);anchor_posi(1,6)-x(1) anchor_posi(2,6)-x(2) Cp^2*t0_s-2*Cp^2*t6*t0+Cp^2*t6^2+f]>=0;
    [1 0 anchor_posi(1,7)-x(1);0 1 anchor_posi(2,7)-x(2);anchor_posi(1,7)-x(1) anchor_posi(2,7)-x(2) Cp^2*t0_s-2*Cp^2*t7*t0+Cp^2*t7^2+f]>=0;
    [1 0 anchor_posi(1,8)-x(1);0 1 anchor_posi(2,8)-x(2);anchor_posi(1,8)-x(1) anchor_posi(2,8)-x(2) Cp^2*t0_s-2*Cp^2*t8*t0+Cp^2*t8^2+f]>=0;
    f>=0;
    [1 0 x(1);0 1 x(2);x(1) x(2) z]>=0;

%     c(1)>=0; %weak constraint
%     c(2)>=0;
%     c(3)>=0;
%     c(4)>=0;
%     c(5)>=0;
%     c(6)>=0;
%     c(7)>=0;
%     c(8)>=0;

    
    t1>=t0;t2>=t0;t3>=t0;t4>=t0;t5>=t0;t6>=t0;t7>=t0;t8>=t0;%really important constrain
    [t0_s t0;t0 1]>=0; %


    (t2-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,1));
    (t3-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,1));
    (t4-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,1));
    (t5-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,1));
    (t6-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,1));
    (t7-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,1));
    (t8-t0)*Cp+(t1-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,1));

    (t1-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,2));
    (t3-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,2));
    (t4-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,2));
    (t5-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,2));
    (t6-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,2));
    (t7-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,2));
    (t8-t0)*Cp+(t2-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,2));

    (t2-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,3));
    (t1-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,3));
    (t4-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,3));
    (t5-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,3));
    (t6-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,3));
    (t7-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,3));
    (t8-t0)*Cp+(t3-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,3));

    (t2-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,4));
    (t3-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,4));
    (t1-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,4));
    (t5-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,4));
    (t6-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,4));
    (t7-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,4));
    (t8-t0)*Cp+(t4-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,4));

    (t2-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,5));
    (t3-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,5));
    (t4-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,5));
    (t1-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,5));
    (t6-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,5));
    (t7-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,5));
    (t8-t0)*Cp+(t5-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,5));

    (t2-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,6));
    (t3-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,6));
    (t4-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,6));
    (t5-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,6));
    (t1-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,6));
    (t7-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,6));
    (t8-t0)*Cp+(t6-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,6));

    (t2-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,7));
    (t3-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,7));
    (t4-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,7));
    (t5-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,7));
    (t6-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,7));
    (t1-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,7));
    (t8-t0)*Cp+(t7-t0)*Cp>=norm(anchor_posi(:,8)-anchor_posi(:,7));

    (t2-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,2)-anchor_posi(:,8));
    (t3-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,3)-anchor_posi(:,8));
    (t4-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,4)-anchor_posi(:,8));
    (t5-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,5)-anchor_posi(:,8));
    (t6-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,6)-anchor_posi(:,8));
    (t7-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,7)-anchor_posi(:,8));
    (t1-t0)*Cp+(t8-t0)*Cp>=norm(anchor_posi(:,1)-anchor_posi(:,8));

    cvx_end;
    %%%%%%%%%%%%%%%%%%%%%%%%%% end cvx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%