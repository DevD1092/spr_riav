# spr_riav
As researchers, we firmly believe in the concept of reproducible research.

Thus, we release the codes of simulation, algorithm development, and testing online on GitHub for the research paper "A rotation-invariant additive vector sequence based star pattern recognition" (https://ieeexplore.ieee.org/abstract/document/8430566/). This repository contains the codes for simulation and testing of the star identification algorithm based on a rotation-invariant additive vector sequence.

Please cite our paper (citation details below) properly, if you use the code for simulation, development or testing of a new star identification algorithm or comparison & benchmarking with the above-mentioned alogirhtm.

This code is only for academic and research purposes. Commerical use of this code is not permitted.

Citation of the paper "A rotation-invariant additive vector sequence based star pattern recognition":
@article{samirbhai2018rotation,<br />
  title={A rotation-invariant additive vector sequence based star pattern recognition},<br />
  author={Samirbhai, Mehta Deval and Chen, Shoushun and Low, Kay Soon},<br />
  journal={IEEE Transactions on Aerospace and Electronic Systems},<br />
  year={2018},<br />
  publisher={IEEE}<br />
}


## Repository details

This repository contains four folders. The information about each of the folder is described as below.

#### 1. simulate - Codes for simulating the star images
&nbsp;&nbsp;&nbsp;&nbsp; Convert_Axis_2_AttitudeMatrix.m -- For converting the ECI (Star position in the catalog) frame to the camera frame (Star sesnsor).<br />
&nbsp;&nbsp;&nbsp;&nbsp;Find_neighbor_star_FOV.m -- For finding the number and position of the neighboring stars in a specified FOV from the center star.<br />
&nbsp;&nbsp;&nbsp;&nbsp;PSF.m -- Point Spread Function simulation of the star amongst the pixels.<br />
&nbsp;&nbsp;&nbsp;&nbsp;Plot_sky_images.m -- For simulating the star images at a specific RA & DEC angle along with a defined FOV (this function is used by the Testing technique eventually).<br />
&nbsp;&nbsp;&nbsp;&nbsp;centroider.m -- Finding the centroid of the stars in the image.<br />
&nbsp;&nbsp;&nbsp;&nbsp;SKY2000_Magnitude6_doublestars_0.12.txt -- Star catalog (adopted from SAO) containing the star ID and it's corresponding RA, DEC and Mv information. Stars having a relative magnitude threshold (Mv) of less than 6.0 are selected for making this catalog.<br />
  
#### 2. Testing - Testing as well as implementation code of the proposed star identification algorithm
&nbsp;&nbsp;&nbsp;&nbsp;Testing_proposed_technique.m -- Code for testing and implementation of the proposed technique in all the sceanrios(i.e. positional deviation, false stars or magnitude uncertainty). The scenario specific parameters can be changed in the script.<br />
###### NOTE: The testing and implementation scripts utilize the simulation scripts as well as the input from the SPD directories. So, please change the path of this input accordingly.

#### 3. SPD - Generating the SPD for the propsoed technique
&nbsp;&nbsp;&nbsp;&nbsp;SPD_generate.m -- For generating the SPD for the propsoed technique. Specific parameters (such as the FOV, pixel size, Mv, bin_size, etc.) can be changed inside the script.<br />
&nbsp;&nbsp;&nbsp;&nbsp;SPD_top_4_dist_list_Mv_6.txt -- SPD containing the distances of the nearest four stars to the reference star.<br />\
&nbsp;&nbsp;&nbsp;&nbsp;SPD_vect_patt_Mv_6_dist_1.txt -- SPD of the rotation-invariant additive vector sequence considering the nearest star as the starting point of the vector pattern.<br />
&nbsp;&nbsp;&nbsp;&nbsp;SPD_vect_patt_Mv_6_dist_2.txt -- SPD of the rotation-invariant additive vector sequence considering the second nearest star as the starting point of the vector pattern.<br />
&nbsp;&nbsp;&nbsp;&nbsp;SPD_vect_patt_Mv_6_dist_3.txt -- SPD of the rotation-invariant additive vector sequence considering the third nearest star as the starting point of the vector pattern.<br />
&nbsp;&nbsp;&nbsp;&nbsp;SPD_vect_patt_Mv_6_dist_4.txt -- SPD of the rotation-invariant additive vector sequence considering the fourth nearest star as the starting point of the vector pattern.<br />
