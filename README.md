This repository contains the scripts and functions of the TurbIFA model, previously 
published IFA data (White et al., 2018; Thirumalai et al., 2019; Stainbank et al., 2021) 
used as case studies and the code routines used for the sensitivity experiments (Figures 2-4). 
All codes have been generated in MATLAB. The manuscript presenting TurbIFA has been submitted
to Paleoceanography and Paleoclimatology.

IMPORTANT! Monthly surface air temperature and decadal sea surface salinity TraCE-21ka data
(Trace_SAT.mat, SSS_Trace.mat) and the three videotutorials are not included in this repository 
(too heavy). You can download them at TurbIFA Zenodo repository: https://doi.org/10.5281/zenodo.110912861.

1. TurbIFA model
* TurbIFA_Tv3_1.m: MATLAB code that performs a bioturbation simulation to calculate the uncertainty
in Mg/Ca IFA datasets sampled in a bioturbated sediment archive;
* TurbIFA_d18Ov3_1.m: MATLAB code that performs a bioturbation simulation to calculate uncertainty
in δ18O IFA datasets sampled in a bioturbated sediment archive;
* quantifaerrv3.m: MATLAB function that calculates the uncertainty in Mg/Ca and δ18O IFA datasets
sampled in a non-bioturbated and bioturbated sediment archive;
* dynheatmap.m MATLAB function that performs a bioturbation simulation using as input a longer
 age time series to avoid biases due to e.g., high SMLD and low SAR.
* analiBH.m: MATLAB function that calculate the probability of picking individuals that grew during
time slices older than mean age of the IFA time slice. It samples “foraminiferal individuals”
from an exponential age distribution using SMLD/SAR as an scale parameter
* Trace_SAT.mat: Global monthly resolved, surface air temperature (SAT) time series spanning
the last 21,000 years (TraCE-21ka, He, 2011);
* SSS_Trace.mat: Global decadally resolved, sea surface salinity (SSS) time series spanning
the last 21,000 years (TraCE-21ka, He, 2011);
* Trace_coordinates.mat: Latitudinal and longitudinal values of SAT and SSS TraCE-21ka outputs (He, 2011);
* datafor_d18Osw.mat: Contains sea-level reconstruction (Lambeck et al.,2014) and the δ18Osw-S
slopes and intercepts of Legrande & Schmidt (2006).

2. Other scripts
*   bioprobv3_1.m: Matlab code to calculate the probability of picking individuals that grew
   during time slices older than mean age of the IFA time slice in e.g., datasets older than
   Trace-21ka outputs (LGM).

3. Sensitivity experiments
* Figure2_code.m: Matlab code used to generate manuscript Figure 2.
* Figure3_code.m: Matlab code used to generate manuscript Figure 3.
* Figure4_code.m: Matlab code used to generate manuscript Figure 4.
* Fig2_4_Tseries.mat: Synthetic SST time series used to generate manuscript Figures 2-4

4. IFA datasets for manuscript case examples (IFAdata_casexamples.mat)
* LH_data1.mat: Mg/Ca (ºC) IFA data, core top (example 1, core MGL1208-14MC/12GC: Globigerinoides ruber);
* Past_data1.mat: Mg/Ca (ºC) IFA data, Younger Dryas (example 1, core MGL1208-14MC/12GC, G. ruber);
* LH_data2.mat: δ18O (‰) IFA data, core top (example 2, core SO189-119KL, G. ruber);
* Past_data2.mat: δ18O (‰) pseudo-IFA data, early Holocene (example 2, core SO189-119KL, G. ruber);
* LH_data3.mat: δ18O (‰) IFA data, core top (example 3, IODP Site U1467A, G. ruber)
* Past_data3.mat: δ18O (‰) IFA data, MIS 9e (example 3, IODP Site U1467C, G. ruber).

5.TurbIFA examples
* Example1_TurbIFA_Tv3_1.m
* Example2_TurbIFA_d18Ov3_1.m
* data_examplesTurBIFA.mat

6. User Guide
* TurbIFA_User_Guide.pdf

7. Video-tutorials
* videotutorial_TurbIFA_T_constant.mp4
* videotutorial_TurbIFA_T_dynamical.mp4
* videotutorial_bioprob.mp4

