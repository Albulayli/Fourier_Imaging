# Fourier_Imaging

**** Using Our Plane-Wave Stolt's and Slant-Stack Migration Methods with PICMUS (IUS-2016 Version) ****
Reference: M. Albulayli and D. Rakhmatov, "Fourier-Domain Depth Migration for Plane-Wave Ultrasound Imaging",
IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control (DOI:10.1109/TUFFC.2018.2837000).


1. File ReferenceDAS.zip should be unzipped first, then file pwFDM.zip should be unzipped into the same location.
(NOTE: See "ReferenceDAS-README.txt" for PICMUS-related details.)


2. File pwFDM.zip contains "pwFDM" folder, which includes the following image reconstruction scripts:

[2a] "script_pwStolt_expcon.m" - application of our plane-wave Stolt's migration method to the "expcon.mat" and
"expcon_scan.mat" inputs. 
(NOTES: It creates "rf_scan" and "dataset" objects from the input MAT files "../Reference/expcon_scan.mat" and
"../Reference/expcon.mat", respectively.  It also creates the appropriate "settings" structure [lines 58-64,72]
and forward-shifts "RFdata = dataset.data(:,:,pw)" along the time-axis by "lense_shift = 32" samples [lines 71,73, 
75].  Individual beamformed data frames "BDcomp" are produced by calling the function "pwStolt(RFdata,settings)"
[line 77], which is followed by backward-shifting that data along the z-axis by "lense_shift" samples [line 78].
Image evaluation is performed by calling the function "exec_evaluation_contrast_speckle_exp('expcon_image_pwStolt.mat',
'expcon_result_pwStolt.txt')" that can be found in the "Reference" folder.  The MAT file "expcon_image_pwStolt.mat"
stores the "expcon_image" structure [lines 102-107] with "expcon_image.scan = scan" and "expcon_image.data = 
resampled_envelope_beamformed_data".  Our "scan" structure [lines 36-53] has the same x-axis and z-axis as those used
during image reconstruction, except that our z-axis is limited by "min(rf_scan.z_axis)" and "max(rf_scan.z_axis)".
Our "scan.x_axis" and "scan.z_axis" are used to obtain "resampled_envelope_beamformed_data" [lines 89-96].  2D linear
interpolation of "envelope_beamformed_data" can be controlled by changing "scan.x_axis" and "scan.z_axis" [lines 36-37,
43-44]; in this version, the script simply crops full-sized 3328x128 envelope data to 1216x128 in size, covering 5-50
mm of imaging depth.)

[2b] "script_pwStolt_expres.m" - application of our plane-wave Stolt's migration method to the "expres.mat" and
"expres_scan.mat" inputs.
(NOTES: See [2a] above, replacing "expcon" with "expres".  Image evaluation is performed by calling the function 
"exec_evaluation_resolution_distorsion_exp('expres_image_pwStolt.mat','expres_result_pwStolt.txt')" [line 168].)

[2c] "script_pwStolt_vcross.m" - application of our plane-wave Stolt's migration method to the "vcross.mat" and
"vcross_scan.mat" inputs.
(NOTES: See [2a] above, replacing "expcon" with "vcross".  The full-sized RF data frames are 1536x128.)

[2d] "script_pwStolt_vlong.m" - application of our plane-wave Stolt's method to the "vlong.mat" and "vlong_scan.mat"
inputs.
(NOTES: See [2a] above, replacing "expcon" with "vlong".  The full-sized RF data frames are 1536x128.)

[2e] "script_pwSlantStack_expcon.m" - application of our plane-wave slant-stack migration method to the "expcon.mat"
and "expcon_scan.mat" inputs. 
(NOTES: See [2a] above, replacing "pwStolt" with "pwSlantStack".  Additional "settings" parameters for setting up
slant samples and their spacing are defined on lines 65-67.  Image evaluation is performed by calling the function
"exec_evaluation_contrast_speckle_exp('expcon_image_pwSlantStack.mat','expcon_result_pwSlantStack.txt')" [line 171].)

[2f] "script_pwSlantStack_expres.m" - application of our plane-wave slant-stack migration method to the "expres.mat"
and "expres_scan.mat" inputs.
(NOTES: See [2a] above, replacing "expcon" with "expres" and replacing "pwStolt" with "pwSlantStack".  Additional
"settings" parameters for setting up slant samples and their spacing are defined on lines 65-67.  Image evaluation
is performed by calling the function "exec_evaluation_resolution_distorsion_exp('expres_image_pwSlantStack.mat',
'expres_result_pwSlantStack.txt')" [line 171].)

[2g] "script_pwSlantStack_vcross.m" - application of our plane-wave slant-stack migration method to the "vcross.mat"
and "vcross_scan.mat" inputs.
(NOTES: See [2a] above, replacing "expcon" with "vcross" and replacing "pwStolt" with "pwSlantStack".  Additional
"settings" parameters for setting up slant samples and their spacing are defined on lines 65-67.  The full-sized
RF data frames are 1536x128.)

[2h] "script_pwSlantStack_vlong.m" - application of our plane-wave slant-stack migration method to the "vlong.mat"
and "vlong_scan.mat" inputs.
(NOTES: See [2a] above, replacing "expcon" with "vlong" and replacing "pwStolt" with "pwSlantStack".  Additional
"settings" parameters for setting up slant samples and their spacing are defined on lines 65-67.  The full-sized
RF data frames are 1536x128.)


3. Acknowledgements:
We sincerely thank PICMUS developers (https://www.creatis.insa-lyon.fr/Challenge/IEEE_IUS_2016/organizers) for making
their experimental data and MATLAB code public.  We would greatly appreciate any feedback from potential users.
