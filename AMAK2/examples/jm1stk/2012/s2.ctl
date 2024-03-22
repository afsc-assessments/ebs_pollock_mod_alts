# jjm.dat file (controls actual data and configurations)   
south.dat   
South_and_offshore   
# Selectivity sharing vector (number_fisheries + number_surveys)   
#Fsh 1 Fsh_3 Fsh_4 Srv_1 Srv_2 Srv_3 Srv_6 Srv_7 Srv_8  
#N_Chile_Fshry CS_Chile_Fishery International Chile_AcousCS Chile_AcousN Chilean_CPUE DEPM Chinese_CPUE EU_CPUE USRR  
# 2 3 1 2 3 4 5 6 7  
1 1 1 2 1 1 2 1 1 1  
1 2 3 1 1 2 4 3 3 3  
#Sr_type 
2   
#AgeError 
0   
#Retro 
0   
#Steepness 
0.8   0.2   -6   
#SigmaR 
0.6   0.2   7   
#yrs_sr 
1977   2011   
#Linf
76.464 0.1 -4   
#K
0.09 0.1 -4   
#Lo_Len
18 0.1 -4   
#Sigma_len
0.09 0.1 -4   
#Natural_Mortality 
0.23   0.1   -8   
# NEW npars_mage
0
# NEW Mage_in

# phase_Mage
-5
#Phase_Random_walk_M 
-4
#Nyrs_Random_walk_M 
0
#Random_walk_M_yrs blank if nyrs==0

#Random_walk_M_sigmas blank if nyrs==0

#catchability 
0.5000   0.0200   0.0600   0.3500   0.0600   0.0600   0.0600   
1.2   1.2   1.2   1.2   1.2   1.2   1.2   
3   2   3   3   1   1   1   
#q_power                    
1   1   1   1   1   1   1   
1.2   1.2   1.2   1.2   1.2   1.2   1.2   
-1   -1   -1   -1   -1   -1   -1   
#Random_walk_q_phases                    
-1  -1  -1  -1  -1  -1  -1  
#Nyrs_Random_walk_q
0  0  0  0  0  0  0  
#Random_walk_q_yrs blank if nyrs==0

#Random_walk_q_sigmas blank if nyrs==0

#q_agemin                    
2   2   2   2   2   2   2   
#q_agemax                    
10   10   10   10   10   10   10   
#junk                    
0.05   
#n_proj_yrs                    
10   
#---------------------------------------------------------
# Fishery 1 N Chile  
1  #selectivity type
9  #n_sel_ages
2  #phase sel
1  #curvature penalty
25  #Dome-shape penalty
# nYears of selectivity change Fishery 1 N Chile  
3
# Years of selectivity change Fishery 1 N Chile  
1987 1990 2003
0.7 0.7 0.7
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Fishery 2, Central South Chile  
1  
10  
3  
1  
12.5  
# nYears of selectivity change Fishery 2, Central South Chile  
5 
# Years of selectivity change Fishery 2, Central South Chile  
1987 1992 1996 2000 2004
0.7  0.7  0.7 0.7 0.7
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Fishery 3 International  
1  
9  
3  
1  
12.5  
# NYears of selectivity change Fishery 3 International  
4
#Actual years
1996 2000 2005 2011
0.7  0.7 0.7 0.7
# Initial values for coefficitients at each change (one for every change plus 1)  
# 2 3 4 5 6 7 8 9 10 11 12  
0.2 0.7 1 1 1 1 1 1 1 1 1 1
#---------------------------------------------------------  
# Index number 1 AcousCS  
1  
10  
2  
0.25  
100  
1
2005
0.7
# Initial values for coefficitients at each change (one for every change plus 1)  
0.3 1 1 1 1 1 1 1 1 1 1 1
#---------------------------------------------------------  
# Index number 2 Acous_N  
1  
10  
-5  
0.25  
100  
0 
0.3 1 1 1 1 1 1 1 1 1 1 1 
#--------------------------------------------------------- 
# Index number 3 Chile_CPUE  
1  
10  
-5  
0.25  
100  
0 
0.8 1 1 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 4 DEPM  
1  
10  
3  
0.25  
100 
1
2003
0.7
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 5 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 6 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
# Index number 7 EU_CPUE  
1  
10  
-5  
0.25  
100  
0 #0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.04 0.43 0.93 1 1 1 1 1 1 1 1 1  
#---------------------------------------------------------
#Test  
123456789  
