This repository saves the source codes used for generating the data in a paper 
"Molecular clockwork hypothesis for the KaiABC circadian oscillations" 
by Masaki Sasai and Shin Fujishiro

[I] Data for Figure 2 representing the efffect of KaiA sequestration:
(1) Amplitude_Sequest.f90: The phase diagram showing the amplitude of oscillations was calculated.
(2) Trajectories_Sequest.f90: Example trajectories were calculated.
 
[II] Data for Figure 3 comparing the numerical results obtained by the Langevin equation and the mean-field theory.
(1) critline.f90: Identifying the crtical line in the mean-field theory.
(2) tamp2.f90: The pahse diagram plotted on the plane varying coupling strenghts was obtained by calculating the amplitude of oscillations.
(3) s3.f90:Amplitudes, variances, and distributions were calculated using both the Langevin equation and the mean-field theory.

[III] Data for Figure 4
(1) s3.f90:Oscillation frequency was calculated from the FFT results of the calculated trajectories.