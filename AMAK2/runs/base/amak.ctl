#dataFile
pm_2023.dat
#modelName
ebswp
#nStocks
1
#nameStocks
EBS_Pollock
#selectivity_shar_matrix 
1 1 1 1 1 1
1 2 2 2 2 2 
1 1 2 3 4 2
#nregbyStock
1
#Sr_type 
1
#AgeError 
0
#Retro 
0
#recMatrix
1
#Steepness 
0.7 
0.1  
-6
#SigmaR 
0.67 
0.05 
-4
#phase_Rzero 
4 
#Nyrs_sr 
44
#yrs_sr_1 
1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993
1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009
2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021
#GrowMatrix 
1 
#Linf
74.4	
0.1	
-4																																									
#K
0.16	
0.1	
-4																																									
#Lo_Len
18	
0.1	
-4																																									
#Sigma_len
0.09 
0.1	
-4																																									
#NMatrix
1
#Natural_Mortality 
0.3 
0.2 
-4
# npars_mage
3
# ages_M_changes
1 2 3
# NEW Mage_in
0.9 0.6 0.3
#phase_Mage
-5
#phase_rw_Mage
-5
#Nyrs_Random_walk_M 
0
#Random_walk_M_yrs blank if nyrs==0

#Random_walk_M_sigmas blank if nyrs==0

#catchability 
1 1 1 1 .1
#catchability prior sigma
.8 .8 10 10 10
#q_phases 
4 4 1 1 3
#q_power 
1 1 1 1 1
#q_power_CV
.1 .1 .1 .1 .1
#q_power_phase
-1 -1 -1 -1 -1
#Random_walk_q_phases 
-4 -4 -4 -4 -4
#Nyrs_Random_walk_q
0 0 0 0 0
#Random_walk_q_yrs blank if nyrs==0

#Random_walk_q_sigmas blank if nyrs==0

#q_agemin 
5 2 1 1  2
#q_agemax 
12 5 1 1 5
#use_vb_wt_age 
0
#n_proj_yrs 
20
#Fsh_selopt_1 
1
#Fsh_nages_1 
9
#Fsh_ph_1 
3
#Fsh_curvpen_1 
1
#Fsh_domepen_1 
1
#Fsh_sel_change_1 
# nyrs fish selectivity changes
58
# Yrs fish selectivity changes
1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981
1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998
1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015
2016 2017 2018 2019 2020 2021 2022
# sigm fish selectivity changes
0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 
# Initial values for coefficitients at each change (one for every change plus 1) 
0 0.1 0.5 0.8 1 1 1 1 1 1 1 1 1 1 1
#
#BTS_selopt 
1
#Ind_nages
8
#Ind_ph
3
#Ind_curvPen
0.1
#Ind_domePen
0.1
# nyrs fish selectivity changes
40
# Yrs fish selectivity changes
1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998
1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015
2016 2017 2018 2019 2020 2021 2022
# sigm fish selectivity changes
0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
# Ages of selectivity change 
0 0.1 0.5 0.8 1 1 1 1 1 1 1 1 1 1 1
#ATS_selopt 
1
#Ind_nages
8
#Ind_ph
-3
#Ind_curvPen
0.1
#Ind_domePen
0.1
# nyrs fish selectivity changes
0
# Yrs fish selectivity changes

# sigm fish selectivity changes

# Ages of selectivity change ATS
0 0.9 0.8 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7
#Age_1_BTS_selopt 
1
#Ind_nages
8
#Ind_ph
-3
#Ind_curvPen
0.1
#Ind_domePen
0.1
# nyrs fish selectivity changes
0
# Yrs fish selectivity changes

# sigm fish selectivity changes

# Ages of selectivity change Fishery 2
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
#Age_1_ATS_selopt 
1
#Ind_nages
8
#Ind_ph
-3
#Ind_curvPen
0.1
#Ind_domePen
0.1
# nyrs fish selectivity changes
0
# Yrs fish selectivity changes

# sigm fish selectivity changes

# Ages of selectivity change Fishery 2
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
#AVO_selopt 
1
#Ind_nages
8
#Ind_ph
-3
#Ind_curvPen
0.1
#Ind_domePen
0.1
# nyrs fish selectivity changes
0
# Yrs fish selectivity changes

# sigm fish selectivity changes

# Ages of selectivity change Fishery 2
0 0.9 0.8 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7
#Age_1_BTS_selopt 
#pop_wtage																																																											
0.084881666	0.195868126	0.3542318	0.47891	0.5902872	0.698033	0.8060948	0.9029996	0.9874292	1.0689416	1.1511662	1.2370584	1.3379166	1.4536958	1.5752528																																													
#mat_age																																																											
0 0.008 0.289 0.641 0.842 0.901 0.947 0.963 0.97 1 1 1 1 1 1                                        																																																											
#Test
123456789
