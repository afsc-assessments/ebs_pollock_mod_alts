# jjm.dat file (controls actual data and configurations)  
North2.dat  
N3  
# Selectivity sharing vector (number_fisheries + number_surveys)  
#Fsh 1 Srv_1 Srv_2 
#Peruvian Acoustic_Peru Peru_CPUE #  
# 2 3 4 1 2 3 4 5 6 7 8 9  
1 2 1
1 1 1
#Sr_type 
2  
#AgeError 
0  
#Retro 
0  
#Steepness 
0.8 300 -6  
#SigmaR 
0.6 15 -4  
#yrs_sr 
1977 2010 
#Linf
80.4	0.1	-4																																									
#K
0.16	0.1	-4																																									
#Lo_Len
18.	0.1	-4																																									
#Sigma_len
0.09 0.1	-4																																									
#Natural_Mortality 
0.33 0.05 -4  
# NEW npars_mage
0
# NEW Mage_in

# phase_Mage
-5
#Phase_Random_walk_M 
4
#Nyrs_Random_walk_M 
1
#Random_walk_M_yrs blank if nyrs==0
2000
#Random_walk_M_sigmas blank if nyrs==0
0.3
#catchability 
1.0  0.5  
0.15  1.2  
-3  4  
#q_power                    
1  1  
1.2  1.2  
-1  -1  
#Random_walk_q_phases                    
-3  -1 
#Nyrs_Random_walk_q
0  0  
#Random_walk_q_yrs blank if nyrs==0
# 2000
#Random_walk_q_sigmas blank if nyrs==0
# 10.
#q_agemin                    
1  1 
#q_agemax                    
6  6 
#junk                    
0.05  
#n_proj_yrs                    
10  
#---------------------------------------------------------
# Fishery 1 Peru  
1  
10  
4  
0.7  
12.5  
# Years of selectivity change Fishery 3 Peru  
1 
2002
0.7
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 .1 .1 .1 .1 .1 .1 .01 .01 .01  
#---------------------------------------------------------
# Index number 2 Acoustic_Peru  
1  
10  
-5  
0.25  
100  
0 
0.8 0.93 1.0 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 6 Peru_CPUE  
1  
10  
-5  
0.25  
100  
0 
0.2 0.7 1 .1 .1 .1 .1 .1 .1 .01 .01 .01  
#---------------------------------------------------------
#Test  
123456789  
