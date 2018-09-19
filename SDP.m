function [x]=SDP(d, anchor_posi, c_speed)
% Implementation of paper "Efficient convex relaxation methods for robust 
%    target localization by a sensor network using time difference of
%    arrivals, K. Yang, G. Wang, and Z. Luo, 2009


% coded by Zhenqinag Su 

%% 
d1_observe = d(1);
d2_observe = d(2);
d3_observe = d(3);
d4_observe = d(4);
d5_observe = d(5);
d6_observe = d(6);
d7_observe = d(7);
d8_observe = d(8);



c=c_speed;%c=299792458;%m/s
tao12=d1_observe-d2_observe;
tao13=d1_observe-d3_observe;
tao14=d1_observe-d4_observe;
tao15=d1_observe-d5_observe;
tao16=d1_observe-d6_observe;
tao17=d1_observe-d7_observe;
tao18=d1_observe-d8_observe;
tao23=d2_observe-d3_observe;
tao24=d2_observe-d4_observe;
tao25=d2_observe-d5_observe;
tao26=d2_observe-d6_observe;
tao27=d2_observe-d7_observe;
tao28=d2_observe-d8_observe;
tao34=d3_observe-d4_observe;
tao35=d3_observe-d5_observe;
tao36=d3_observe-d6_observe;
tao37=d3_observe-d7_observe;
tao38=d3_observe-d8_observe;
tao45=d4_observe-d5_observe;
tao46=d4_observe-d6_observe;
tao47=d4_observe-d7_observe;
tao48=d4_observe-d8_observe;
tao56=d5_observe-d6_observe;
tao57=d5_observe-d7_observe;
tao58=d5_observe-d8_observe;
tao67=d6_observe-d7_observe;
tao68=d6_observe-d8_observe;
tao78=d7_observe-d8_observe;

tao_wave=[tao12 tao13 tao14 tao15 tao16 tao17 tao18 tao23 tao24 tao25 tao26 tao27 tao28 tao34 tao35 tao36 tao37 tao38 tao45 tao46 tao47 tao48 tao56 tao57 tao58 tao67 tao68 tao78]';
tao_wave=tao_wave/c;%transfer to time delay.
tao_line=-tao_wave;
G=[1 -1 0 0 0 0 0 0
    1 0 -1 0 0 0 0 0
    1 0 0 -1 0 0 0 0
    1 0 0 0 -1 0 0 0
    1 0 0 0 0 -1 0 0
    1 0 0 0 0 0 -1 0
    1 0 0 0 0 0 0 -1
    0 1 -1 0 0 0 0 0
    0 1 0 -1 0 0 0 0
    0 1 0 0 -1 0 0 0
    0 1 0 0 0 -1 0 0
    0 1 0 0 0 0 -1 0
    0 1 0 0 0 0 0 -1
    0 0 1 -1 0 0 0 0
    0 0 1 0 -1 0 0 0
    0 0 1 0 0 -1 0 0
    0 0 1 0 0 0 -1 0
    0 0 1 0 0 0 0 -1
    0 0 0 1 -1 0 0 0
    0 0 0 1 0 -1 0 0
    0 0 0 1 0 0 -1 0
    0 0 0 1 0 0 0 -1
    0 0 0 0 1 -1 0 0
    0 0 0 0 1 0 -1 0
    0 0 0 0 1 0 0 -1
    0 0 0 0 0 1 -1 0
    0 0 0 0 0 1 0 -1
    0 0 0 0 0 0 1 -1
    ];
F=[2*G'*G G'*(tao_line-tao_wave);(tao_line-tao_wave)'*G tao_line'*tao_line+tao_wave'*tao_wave];
delta=0.000001;

%%%%%%%%%%%%%%%%%%%%%%%%%% cvx %%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_begin sdp quiet
%cvx_solver sedumi
variable x(2,1)
variable t(8,1)
variable T(8,8)
variable z

nn=[T, t;t', 1];

minimize trace(nn*F) + delta*(sum(sum(T)))


subject to

T(1,1)==c^(-2)*[anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,1) anchor_posi(2,1) -1]');
T(2,2)==c^(-2)*[anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,2) anchor_posi(2,2) -1]');
T(3,3)==c^(-2)*[anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,3) anchor_posi(2,3) -1]');
T(4,4)==c^(-2)*[anchor_posi(1,4) anchor_posi(2,4) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,4) anchor_posi(2,4) -1]');
T(5,5)==c^(-2)*[anchor_posi(1,5) anchor_posi(2,5) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,5) anchor_posi(2,5) -1]');
T(6,6)==c^(-2)*[anchor_posi(1,6) anchor_posi(2,6) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]');
T(7,7)==c^(-2)*[anchor_posi(1,7) anchor_posi(2,7) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]');
T(8,8)==c^(-2)*[anchor_posi(1,8) anchor_posi(2,8) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]');

T(1,2)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,2) anchor_posi(2,2) -1]'));
T(1,3)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,3) anchor_posi(2,3) -1]'));
T(1,4)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,4) anchor_posi(2,4) -1]'));
T(1,5)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,5) anchor_posi(2,5) -1]'));
T(1,6)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]'));
T(1,7)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(1,8)>=c^(-2)*abs([anchor_posi(1,1) anchor_posi(2,1) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));

T(2,3)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,3) anchor_posi(2,3) -1]'));
T(2,4)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,4) anchor_posi(2,4) -1]'));
T(2,5)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,5) anchor_posi(2,5) -1]'));
T(2,6)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]'));
T(2,7)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(2,8)>=c^(-2)*abs([anchor_posi(1,2) anchor_posi(2,2) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));

T(3,4)>=c^(-2)*abs([anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,4) anchor_posi(2,4) -1]'));
T(3,5)>=c^(-2)*abs([anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,5) anchor_posi(2,5) -1]'));
T(3,6)>=c^(-2)*abs([anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]'));
T(3,7)>=c^(-2)*abs([anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(3,8)>=c^(-2)*abs([anchor_posi(1,3) anchor_posi(2,3) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));

T(4,5)>=c^(-2)*abs([anchor_posi(1,4) anchor_posi(2,4) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,5) anchor_posi(2,5) -1]'));
T(4,6)>=c^(-2)*abs([anchor_posi(1,4) anchor_posi(2,4) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]'));
T(4,7)>=c^(-2)*abs([anchor_posi(1,4) anchor_posi(2,4) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(4,8)>=c^(-2)*abs([anchor_posi(1,4) anchor_posi(2,4) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));

T(5,6)>=c^(-2)*abs([anchor_posi(1,5) anchor_posi(2,5) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,6) anchor_posi(2,6) -1]'));
T(5,7)>=c^(-2)*abs([anchor_posi(1,5) anchor_posi(2,5) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(5,8)>=c^(-2)*abs([anchor_posi(1,5) anchor_posi(2,5) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));

T(6,7)>=c^(-2)*abs([anchor_posi(1,6) anchor_posi(2,6) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,7) anchor_posi(2,7) -1]'));
T(6,8)>=c^(-2)*abs([anchor_posi(1,6) anchor_posi(2,6) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));
T(7,8)>=c^(-2)*abs([anchor_posi(1,7) anchor_posi(2,7) -1]*[1 0 x(1);0 1 x(2);x(1) x(2) z]*([anchor_posi(1,8) anchor_posi(2,8) -1]'));


[T t;t' 1]>=0;
[1 0 x(1);0 1 x(2);x(1) x(2) z]>=0;

cvx_end;
%%%%%%%%%%%%%%%%%%%%%%%%%% end cvx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
