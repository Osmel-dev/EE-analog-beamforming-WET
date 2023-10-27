# EE-analog-beamforming-WET
This compilation of MatLab scripts is used to reproduce the results illustrated in Figs. 2, 4, and 5 in [REF]. The repository contains the following files.

## Figure scripts 
- Fig2.m: allows reproducing Fig.2 "Transfer characteristic of the Powercast P2110B RF-EH circuit"
- Fig4.m: allows reproducing Fig.4 "PB's energy consumption vs maximum charging time"    
- Fig5.m: allows reproducing Fig.5 "PB's energy consumption vs \alpha"

## Auxiliary functions
- allAtOnce.m: This function implements the solution of P3 in Section III.C
- orderAgnostic.m: This function implements the solution of the optimization problem P2 (see Section II.A "Order-agnostic time division")
- orderAware.m: This function implements Algorithm 1 "Order-aware time division".
- channelModelURA.m: This function generates realizations of the wireless channel according to the model in [R1,R2] for a uniform rectangular array
- plotSettings.m: This auxiliary function sets the basic properties of the plots

## References 
- [R1] O. M. Rosabal, et al., _"Energy-Efficient Analog Beamforming for RF-WET with Charging Time Constraint"_, _IEEE Trans. Veh. Technol. (submitted)_, 2023.
- [R2] W. Tan, et al, _"Achievable sum-rate analysis for massive MIMO systems with different array configurations,"_ _WCNC_, New Orleans, LA, USA, 2015, pp. 316-321, doi: 10.1109/WCNC.2015.7127489.
- [R3] C. A. Balanis, _Antenna theory: analysis and design_. John Wiley & Sons, 2015
