# Semidefinite-Programming-SDP-optimization
Codes for paper "Semidefinite Programming for NLOS Error Mitigation in TDOA Localization"

Run the test.m file to simulate the algorithms including SDP, SDP-robust and the proposed SDP. 

The simulation settings here are a bit different from the original paper in two ways (but the algorithm is the same): a. here
the NLOS distributions are randomly generated (how many links and which links are in NLOS); b. The source nodes, whoes position to be 
found, are randomly generated. In the paper, NLOS distributions are emulated and the source positions are evenly distributed but fixed. 

There are several hyper-parameters need to be tuned/optimized to obtain the best performance. However, in the implementation, the hyper-parameters 
are only manually selected and therefore they are not optimized. 





