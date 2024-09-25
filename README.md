# DataDrivenNoisyIO
This is the simulation code for the paper 

*Controller Synthesis from Noisy-Input Noisy-Output Data*

Please see the link 
>https://arxiv.org/abs/2402.02588


## Instructions
*BatchReactorIO.m* is for the first case, Section 5.1

*Case_pl_neq_n.m* is for the second case, Section 5.2

*Comparison_SDP.m* and *Comparison_SOS* is for the third case, Section 5.3



## Additions
You need [CVX](https://cvxr.com/cvx/).

The code of *Comparison_SOS* is from the paper 
[*Superstabilizing Control of Discrete-Time ARX Models under Error in Variables*](https://github.com/Jarmill/eiv_arx).

To run the code *test_Dual_SS.m* in *Comparison_SOS*, you need [YALMIP](https://yalmip.github.io/) and [Mosek](https://www.mosek.com/).
