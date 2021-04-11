# Semidefinite-Programming-SDP-optimization
Codes for paper "Semidefinite Programming for NLOS Error Mitigation in TDOA Localization". 
** For those who don't have access to the paper, the PDF file was also uploaded in this repo with name 'Semidefinite Programming for NLOS Error Mitigation in TDOA Localization.pdf'. 

Run the test.m file to simulate the algorithms including SDP, SDP-robust and the proposed SDP. 

The simulation settings here are a bit different from the original paper in two ways (but the algorithm is the same): a. here
the NLOS distributions are randomly generated (how many links and which links are in NLOS); b. The source nodes, whoes position to be 
found, are randomly generated. In the paper, NLOS distributions are emulated and the source positions are evenly distributed but fixed. 

There are several hyper-parameters need to be tuned/optimized to obtain the best performance. However, in the implementation, the hyper-parameters 
are only manually selected and therefore they are not optimized. 

Simulation with 1000 MC runs:
https://user-images.githubusercontent.com/15931069/45726434-41970800-bb8d-11e8-8f4e-0b0a52acc1f0.png
Note: (the y-label for graphs in the paper and this project should be 'absolute distance error' not MSE)


*If you want to do performance comparison between your algorithm and the one in this paper, please use the file 'SDP_CL_SU.m'; 

**  The algorithm 'sdp_ce.m' is provided by Prof. Gang Wang, from Ningbo University. Please cite his paper if you use this algorithm; 

Contact me at xiamenhai[at]gmail.com if you have questions. 
