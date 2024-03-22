////////////////////////////////////////////////////////////////////
// Naming Conventions:
//
//  GENERAL:
//    styr, endyr begining year and ending year of model (catch data available)
//    nages       number of age groups considered
//    nyrs_        number of observations available to specific data set
//
//  DATA SPECIFIC:
//    catch_bio   Observed catch biomass
//    fsh        fishery data
//
//  Define indices
//    nind        number of indices
//  Index values
//    nyrs_ind      Number of years of index value (annual)
//    yrs_ind        Years of index value (annual)
//    obs_ind        Observed index value (annual)
//    obs_se_ind    Observed index standard errors (annual)
//  Age-comp values
//    nyrs_ind_age  Number of years index age data available
//    yrs_ind_age   Years of index age value (annual)
//    oac_ind       Observed age comp from index
//    n_sample_ind_age    Observed age comp sample sizes from index
//
//    eac_ind       Expected age comp from index
//
//    sel_ind       selectivity for egg production index
//
//    pred_ind ...
//
//    oac_fsh      Observed age comp from index
//    obs_ind_size  Observed size comp from index
//
//    pred_fsh_age    Predicted age comp from fishery data
//    eac_fsh            Expected age comp for fishery data (only years where data available)
//    eac_ ...
//
//    pred_tmp_ind   Predicted index value for trawl index
//
//    sel_fsh    selectivity for fishery                
//  
//     sel_ch indicates time-varying selectivity change
//  
//    Add bit for historical F
//    Added length part for selectivity
//
//////////////////////////////////////////////////////////////////////////////
 // To ADD/FIX:
 //   parameterization of steepness to work the same (wrt prior) for ricker and bholt
 //   splines for selectivity
 //   two projection outputs need consolidation
//////////////////////////////////////////////////////////////////////////////

DATA_SECTION
  !!version_info+="SPRFMO_Jack_Mackerel;May_2018";
  int iseed 
  !! iseed=1313;
  int cmp_no // candidate management procedure
  int nnodes_tmp;
  !!CLASS ofstream mceval("mceval.dat")
  !!long int lseed=iseed;
  !!CLASS random_number_generator rng(iseed);
  
  int oper_mod
  int mcmcmode
  int mcflag

  !! oper_mod = 0;
  !! mcmcmode = 0;
  !! mcflag   = 1;
 LOCAL_CALCS
  write_input_log<<version_info<<endl;
  tmpstring=adprogram_name + adstring(".dat");
  int on=0;
  if ( (on=option_match(argc,argv,"-ind"))>-1)
  {
    if (on>argc-2 | argv[on+1][0] == '-') 
    { 
      cerr << "Invalid input data command line option"
         " -- ignored" << endl;  
    }
    else
    {
      cntrlfile_name = adstring(argv[on+1]);
    }
  }
  else
  {
      cntrlfile_name =   tmpstring;
  }
  if ( (on=option_match(argc,argv,"-om"))>-1)
  {
    oper_mod  = 1;
    cmp_no = atoi(argv[on+1]);
    cout<<"Got to operating model option "<<oper_mod<<endl;
  }
  if ( (on=option_match(argc,argv,"-mcmc"))>-1)
  {
    mcmcmode = 1;
  }
  global_datafile= new cifstream(cntrlfile_name);
  if (!global_datafile)
  {
  }
  else
  {
    if (!(*global_datafile))
    {
      delete global_datafile;
      global_datafile=NULL;
    }
  }
 END_CALCS
 // Read in "name" of this model...
  !! *(ad_comm::global_datafile) >>  datafile_name; // First line is datafile (not used by this executable)
  !! *(ad_comm::global_datafile) >>  model_name; 
  !! ad_comm::change_datafile_name(datafile_name);
  init_int styr
  init_int endyr
  init_int rec_age
  init_int oldest_age
  !! log_input(styr);
  !! log_input(endyr);
  !! log_input(rec_age);
  !! log_input(oldest_age);
//------------LENGTH INTERVALS
  init_int nlength
  init_vector len_bins(1,nlength)
  !! log_input(nlength);
  !! log_input(len_bins);

  int nages
  !!  nages = oldest_age - rec_age + 1;
  int styr_rec
  int styr_sp
  int endyr_sp
  int nyrs
  !! nyrs          = endyr - styr + 1;
  int mc_count;
  !!  mc_count=0;
  !! styr_rec = (styr - nages) + 1;     // First year of recruitment
  !! styr_sp  = styr_rec - rec_age - 1 ;    // First year of spawning biomass  
  vector yy(styr,endyr);
  !! yy.fill_seqadd(styr,1) ;
  vector aa(1,nages);
  !! aa.fill_seqadd(rec_age,1) ;
  int junk;
// Fishery specifics
  init_int nfsh                                   //Number of fisheries
  imatrix pfshname(1,nfsh,1,2)
  init_adstring fshnameread;
 LOCAL_CALCS
  for(k=1;k<=nfsh;k++) 
  {
    pfshname(k,1)=1; 
    pfshname(k,2)=1;
  }    // set whole array to equal 1 in case not enough names are read
  adstring_array CRLF;   // blank to terminate lines
  CRLF+="";
  k=1;
  for(i=1;i<=strlen(fshnameread);i++)
  if(adstring(fshnameread(i))==adstring("%")) {
    pfshname(k,2)=i-1; 
    k++;  
    pfshname(k,1)=i+1;
  }
  pfshname(nfsh,2)=strlen(fshnameread);
  for(k=1;k<=nfsh;k++)
  {
    fshname += fshnameread(pfshname(k,1),pfshname(k,2))+CRLF(1);
  }
  log_input(datafile_name);
  log_input(model_name);
  log_input(styr);
  log_input(endyr);
  log_input(rec_age);
  log_input(oldest_age);
  log_input(nfsh);
  log_input(fshname);
 END_CALCS
  init_matrix catch_bio_in(1,nfsh,styr,endyr)
  init_matrix catch_bio_sd_in(1,nfsh,styr,endyr)   // Specify catch-estimation precision
  // !! for (i=1;i<=nfsh;i++) catch_bio(i) += .01; 
  !! log_input(catch_bio_in);
  !! log_input(catch_bio_sd_in);


//  Define fishery age compositions
  init_ivector nyrs_fsh_age(1,nfsh)
  !! log_input(nyrs_fsh_age);
  init_ivector nyrs_fsh_length(1,nfsh)
  !! log_input(nyrs_fsh_length);
  init_imatrix yrs_fsh_age_in(1,nfsh,1,nyrs_fsh_age)
  !! log_input(yrs_fsh_age_in);
  init_imatrix yrs_fsh_length_in(1,nfsh,1,nyrs_fsh_length)
  !! log_input(yrs_fsh_length_in);
  init_matrix n_sample_fsh_age_in(1,nfsh,1,nyrs_fsh_age)    //Years of index index value (annual)
  !! log_input(n_sample_fsh_age_in);
  init_matrix n_sample_fsh_length_in(1,nfsh,1,nyrs_fsh_length)    //Years of index index value (annual)
  !! log_input(n_sample_fsh_length_in);
  init_3darray oac_fsh_in(1,nfsh,1,nyrs_fsh_age,1,nages)
  !! log_input(oac_fsh_in);
  init_3darray olc_fsh_in(1,nfsh,1,nyrs_fsh_length,1,nlength)
  !! log_input(olc_fsh_in);
  init_3darray wt_fsh(1,nfsh,styr,endyr,1,nages)  //values of weights at age
  !! log_input(wt_fsh);

//  Define indices
  init_int nind                                   //number of indices
  !! log_input(nind);
  int nfsh_and_ind
  !! nfsh_and_ind = nfsh+nind;
  vector Francis_wts(1,nfsh_and_ind)
  imatrix pindname(1,nind,1,2)
  init_adstring indnameread;
 LOCAL_CALCS
  for(int k=1;k<=nind;k++) 
  {
    pindname(k,1)=1; 
    pindname(k,2)=1;
  }    // set whole array to equal 1 in case not enough names are read
  int k=1;
  for(i=1;i<=strlen(indnameread);i++)
  if(adstring(indnameread(i))==adstring("%")) {
    pindname(k,2)=i-1; 
    k++;  
    pindname(k,1)=i+1;
  }
  pindname(nind,2)=strlen(indnameread);
  for(k=1;k<=nind;k++)
  {
    indname += indnameread(pindname(k,1),pindname(k,2))+CRLF(1);
  }
  log_input(indname);
 END_CALCS

//  Index values
  init_ivector nyrs_ind(1,nind)                   //Number of years of index value (annual)
  init_imatrix yrs_ind_in(1,nind,1,nyrs_ind)         //Years of index value (annual)
  init_vector mo_ind(1,nind)                      //Month occur 
  init_matrix obs_ind_in(1,nind,1,nyrs_ind)          //values of index value (annual)
  init_matrix obs_se_ind_in(1,nind,1,nyrs_ind)       //values of indices serrs

  vector ind_month_frac(1,nind)
  !! log_input(nyrs_ind);
  !! log_input(yrs_ind_in);
  !! log_input(mo_ind);
  !! ind_month_frac = (mo_ind-1.)/12.;
  !! log_input(obs_ind_in);
  !! log_input(obs_se_ind_in);
  matrix        corr_dev(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
  matrix        corr_eff(1,nfsh,styr,endyr) //Index standard errors (for lognormal)
  matrix         act_eff(1,nfsh,styr,endyr) //Index standard errors (for lognormal)
  vector              ac(1,nind);

  init_ivector nyrs_ind_age(1,nind)               //Number of years of index value (annual)
  !! log_input(nyrs_ind_age);

  init_ivector nyrs_ind_length(1,nind)
  !! log_input(nyrs_ind_length);

  init_imatrix yrs_ind_age_in(1,nind,1,nyrs_ind_age)  //Years of index value (annual)
  !! log_input(yrs_ind_age_in);

  init_imatrix yrs_ind_length_in(1,nind,1,nyrs_ind_length)
  !! log_input(yrs_ind_length_in);

  init_matrix n_sample_ind_age_in(1,nind,1,nyrs_ind_age)         //Years of index value (annual)
  !! log_input(yrs_ind_age_in);

  init_matrix n_sample_ind_length_in(1,nind,1,nyrs_ind_length)         //Years of index lengths (annual)
  !! log_input(n_sample_ind_length_in);

  init_3darray oac_ind_in(1,nind,1,nyrs_ind_age,1,nages);  //values of Index proportions at age
  init_3darray olc_ind_in(1,nind,1,nyrs_ind_length,1,nlength);
  !! log_input(olc_ind_in);

  !! log_input(oac_ind_in);
  init_3darray  wt_ind(1,nind,styr,endyr,1,nages)      //values of Index proportions at age
  !! log_input(wt_ind);

  vector age_vector(1,nages);
  !! for (j=1;j<=nages;j++)
  !!  age_vector(j) = double(j+rec_age-1);

  //Spawning month-----
  init_number spawnmo
  number spmo_frac
  !! spmo_frac = (spawnmo-1)/12.;

  init_matrix age_err(1,nages,1,nages)
  !! log_input(age_err);

  int s // Index for stock
  int r // Index for regime
  int k // Index for fishery or index
  int i // Index for year
  int j // Index for age
 LOCAL_CALCS
  // Rename data file to the control data section... 
  ad_comm::change_datafile_name(cntrlfile_name);
  *(ad_comm::global_datafile) >>  datafile_name; 
  *(ad_comm::global_datafile) >>  model_name; 
  log_input(cntrlfile_name);
 END_CALCS
  // Stock specifics
  init_int nstk                                   //Number of stocks
  imatrix pstkname(1,nstk,1,2)
  init_adstring stknameread;
 LOCAL_CALCS
  for(s=1;s<=nstk;s++) 
  {
    pstkname(s,1)=1; 
    pstkname(s,2)=1;
  }    // set whole array to equal 1 in case not enough names are read
  // adstring_array CRLF;   // blank to terminate lines
  CRLF+="";
  s=1;
  for(i=1;i<=strlen(stknameread);i++)
    if(adstring(stknameread(i))==adstring("%")) {
      pstkname(s,2)=i-1; 
      s++;  
      pstkname(s,1)=i+1;
    }
  pstkname(nstk,2)=strlen(stknameread);
  for(s=1;s<=nstk;s++)
  {
    stkname += stknameread(pstkname(s,1),pstkname(s,2))+CRLF(1);
  }
 END_CALCS
  // Matrix of selectivity mappings--row 1 is index of stock
  //  row 2 is type (1=fishery, 2=index) and row 3 is index within that type
  //  e.g., the following for 2 fisheries and 4 indices means that index 3 uses fishery 1 selectivities,
  //         the other fisheries and indices use their own parameterization
  //  1 1 2 2 1 2 
  //  1 2 1 2 1 4
  init_imatrix sel_map(1,3,1,nfsh_and_ind) 
  // maps fisheries and indices into sequential sel_map for sharing purposes
  !! log_input(datafile_name);
  !! log_input(model_name);
  !! log_input(nstk);
  !! log_input(stkname);
  !! write_input_log<<"# Map shared selectivity: "<< endl;log_input(sel_map);
  !! projfile_name = cntrlfile_name(1,length(cntrlfile_name)-4) + ".prj";

  
  init_ivector    nreg(1,nstk)
  !! log_input(nreg);
  int nregs
  !! nregs = sum(nreg);
  vector cum_regs(1,nstk)
  !! cum_regs.initialize();
  vector R_guess(1,nregs)
  matrix stk_reg_map(1,2,1,nregs)
 LOCAL_CALCS
  //cum_regs(1) = 0;
  for (s=2;s<=nstk;s++)
  {
    cum_regs(s) = sum(nreg(1,s-1));
  }
  int jj=1;
  for (s=1;s<=nstk;s++)
  {
    for (r=1;r<=nreg(s);r++)
    {
      stk_reg_map(1,jj) = s;
      stk_reg_map(2,jj) = r;
      jj++;
    }
  }
 END_CALCS
  init_int    SrType        // 2 Bholt, 1 Ricker
  !! log_input(SrType);
  init_int use_age_err      // nonzero value means use...
  !! log_input(use_age_err);
  init_int retro            // Retro years to peel off (0 means full dataset)
  !! log_input(retro);
  init_imatrix rec_map(1,nstk,1,nreg)
  int nrec
  ivector phase_mean_rec(1,nregs)
  number R_guess_tmp
  int phase_mean_rec_tmp
  int istk_tmp
  int ireg_tmp
  int iregs_tmp
 LOCAL_CALCS
  nrec = max(rec_map);
  phase_mean_rec = 1;
  if ( (on=option_match(argc,argv,"-piner"))>-1)
  {
    if (on>argc-2 | argv[on+1][0] == '-') 
    { 
      cerr << "Invalid piner command line option"
         " -- ignored" << endl;  
    }
    else
    {
      R_guess_tmp = atof(argv[on+1]);
      phase_mean_rec_tmp = -1;
      cout << "Got to piner "<<R_guess_tmp<<endl;
      if (on<=argc-4)
      {
        if (argv[on+2][0] != '-' & argv[on+3][0] != '-')
        {
          istk_tmp = atoi(argv[on+2]);
          ireg_tmp = atoi(argv[on+3]);
        }
        else
        {
          istk_tmp = 1;
          ireg_tmp = 1;
        }
      }
      else
      {
        istk_tmp = 1;
        ireg_tmp = 1;
      }
      R_guess(cum_regs(istk_tmp) + ireg_tmp) = R_guess_tmp;
      phase_mean_rec(cum_regs(istk_tmp) + ireg_tmp) = phase_mean_rec_tmp;
    }
  }
  else
  {
    phase_mean_rec_tmp = 1;
  }
 END_CALCS
  init_vector steepnessprior(1,nrec)
  init_vector cvsteepnessprior(1,nrec)
  init_ivector    phase_srec(1,nrec)

  init_vector sigmarprior(1,nrec)
  vector log_sigmarprior(1,nrec)
  init_vector cvsigmarprior(1,nrec)
  init_ivector    phase_sigmar(1,nrec)
  !! log_input(sigmarprior);
  !! log_input(cvsigmarprior);
  !! log_input(phase_sigmar);
  init_ivector phase_Rzero(1,nregs)
  init_ivector nrecs_est_shift(1,nregs)
  init_imatrix yr_rec_est(1,nregs,1,nrecs_est_shift)
  imatrix styr_rec_est(1,nstk,1,nreg)
  imatrix endyr_rec_est(1,nstk,1,nreg)
 LOCAL_CALCS
  for (s=1;s<=nstk;s++)
  {
    for (r=1;r<=nreg(s);r++)
    {
      styr_rec_est(s,r) = yr_rec_est(cum_regs(s)+r,1);
      endyr_rec_est(s,r) = yr_rec_est(cum_regs(s)+r,nrecs_est_shift(cum_regs(s)+r));
    }
  }
 END_CALCS
  !! log_input(styr_rec_est);
  !! log_input(endyr_rec_est);
  !! log_input(yr_rec_est);
  init_imatrix reg_shift(1,nstk,1,nreg-1)
  imatrix yy_shift_st(1,nstk,1,nreg)
  imatrix yy_shift_end(1,nstk,1,nreg)

//-----GROWTH PARAMETERS--------------------------------------------------
  init_imatrix growth_map(1,nstk,1,nreg)
  int ngrowth
 LOCAL_CALCS
  ngrowth = max(growth_map);
 END_CALCS
  init_vector  Linfprior(1,ngrowth)
  init_vector  cvLinfprior(1,ngrowth)
  init_ivector phase_Linf(1,ngrowth)
  number lw_a 
  number lw_b
  !! lw_a = 0.007778994e-3 ; // LW parameters from Talcahuano sampling, 2008
  !! lw_b = 3.089248476 ; // LW parameters from Talcahuano sampling, 2008
  vector log_Linfprior(1,ngrowth)
  !! log_Linfprior = log(Linfprior);
  !! log_input(Linfprior)
  !! log_input(cvLinfprior)

  init_vector  kprior(1,ngrowth)
  init_vector  cvkprior(1,ngrowth)
  init_ivector phase_k(1,ngrowth)
  vector log_kprior(1,ngrowth)
  !! log_kprior = log(kprior);
  !! log_input(kprior)
  !! log_input(cvkprior)

  init_vector  Loprior(1,ngrowth)
  init_vector  cvLoprior(1,ngrowth)
  init_ivector phase_Lo(1,ngrowth)
  vector log_Loprior(1,ngrowth)
  !! log_Loprior = log(Loprior);
  !! log_input(Loprior)
  !! log_input(cvLoprior)

  init_vector  sdageprior(1,ngrowth)
  init_vector  cvsdageprior(1,ngrowth)
  init_ivector phase_sdage(1,ngrowth)
  vector log_sdageprior(1,ngrowth)
  !! log_sdageprior = log(sdageprior);
  !! log_input(sdageprior)
  !! log_input(cvsdageprior)

//---------------------------------------------------------------------------
  init_imatrix mort_map(1,nstk,1,nreg)
  int nmort
 LOCAL_CALCS
  nmort = max(mort_map);
 END_CALCS
  int npar
 LOCAL_CALCS
  npar = max(max(nrec,ngrowth),nmort);
 END_CALCS
  // Basic M
  init_vector  natmortprior(1,nmort)
  init_vector  cvnatmortprior(1,nmort)
  init_ivector phase_M(1,nmort)
  !! log_input(natmortprior);
  !! log_input(cvnatmortprior);
  !! log_input(phase_M);

  // age-specific M
  init_ivector npars_Mage(1,nmort)
  init_imatrix ages_M_changes(1,nmort,1,npars_Mage)
  init_matrix  Mage_in(1,nmort,1,npars_Mage)
  init_ivector phase_Mage(1,nmort)
  matrix       Mage_offset_in(1,nmort,1,npars_Mage)
  // convert inputs to offsets from prior for initialization purposes
  !! Mage_offset_in.initialize();
 LOCAL_CALCS
  for (r=1;r<=nmort;r++)
  {
    if (npars_Mage(r)>0)
      Mage_offset_in(r) = log(Mage_in(r) / natmortprior(r));
  }
 END_CALCS
  !! log_input(npars_Mage);
  !! log_input(ages_M_changes);
  !! log_input(Mage_in);
  !! log_input(Mage_offset_in);

  // time-varying M
  init_ivector phase_rw_M(1,nstk)
  init_ivector npars_rw_M(1,nstk)
  init_imatrix yrs_rw_M(1,nstk,1,npars_rw_M)
  init_matrix  sigma_rw_M(1,nstk,1,npars_rw_M)
 LOCAL_CALCS
  log_input(phase_rw_M);
  log_input(npars_rw_M);
  log_input(yrs_rw_M);
  log_input(sigma_rw_M);
 END_CALCS

  init_vector qprior(1,nind)      
  vector log_qprior(1,nind)      
  init_vector cvqprior(1,nind)     
  init_ivector phase_q(1,nind)
  !! log_input(qprior);
  !! log_input(cvqprior);
  !! log_input(phase_q);

  init_vector q_power_prior(1,nind)      
  vector log_q_power_prior(1,nind)      
  init_vector cvq_power_prior(1,nind)     
  init_ivector phase_q_power(1,nind)
  // Random walk definition for indices
  init_ivector phase_rw_q(1,nind)
  init_ivector npars_rw_q(1,nind)
  init_imatrix  yrs_rw_q(1,nind,1,npars_rw_q); // Ragged array
  init_matrix sigma_rw_q(1,nind,1,npars_rw_q); // Ragged array
 LOCAL_CALCS
  log_input(phase_rw_q);
  log_input(npars_rw_q);
  log_input(yrs_rw_q);
  log_input(sigma_rw_q);
 END_CALCS

  init_ivector    q_age_min(1,nind)     // Age that q relates to...
  init_ivector    q_age_max(1,nind)     // Age that q relates to...
  !! log_input(q_age_min);
  !! log_input(q_age_max);
  // Need to map to age index range...
  !! for (k=1;k<=nind;k++) {q_age_min(k) =  q_age_min(k) - rec_age + 1; q_age_max(k) = q_age_max(k) - rec_age + 1;}
  !! log_input(q_age_min);
  !! log_input(q_age_max);

  init_int use_vb_wt_age // Flag to use hard-wired wt-age as function of VB params (not inputs)

  number catchbiomass_pen
  number cv_catchbiomass
  !!cv_catchbiomass = 0.05; // 1./(2*0.05*0.05);
  !!catchbiomass_pen= 200.; // 1./(2*0.05*0.05);

  init_int nproj_yrs

  int styr_fut
  int endyr_fut            // LAst year for projections
  int phase_nosr
  number Steepness_UB
  // !! phase_Rzero =  4;
  !! phase_nosr  = -3;
  
  matrix yy_sr(1,nstk,styr_sp,endyr+nproj_yrs);
 LOCAL_CALCS
  yy_sr = 1;
  for (s=1;s<=nstk;s++)
  {
    for (r=2;r<=nreg(s);r++)
    {
      for (i=reg_shift(s,r-1);i<=endyr+nproj_yrs;i++)
      {
        yy_sr(s,i) = r;
      }
    }
  }
 END_CALCS
  
  // Selectivity controls
  // read in options for each fishery
  // Loop over fisheries and indices to read in data (conditional on sel_options)
  ivector   fsh_sel_opt(1,nfsh)
  ivector phase_sel_fsh(1,nfsh)
  vector   curv_pen_fsh(1,nfsh)
  matrix   sel_slp_in_fsh(1,nfsh,1,nyrs)
  matrix   logsel_slp_in_fsh(1,nfsh,1,nyrs)
  matrix   sel_inf_in_fsh(1,nfsh,1,nyrs)
  vector   logsel_slp_in_fshv(1,nfsh)
  vector   sel_inf_in_fshv(1,nfsh)
  vector   logsel_dslp_in_fshv(1,nfsh)
  vector   sel_dinf_in_fshv(1,nfsh)
  matrix   sel_dslp_in_fsh(1,nfsh,1,nyrs)
  matrix   logsel_dslp_in_fsh(1,nfsh,1,nyrs)
  matrix   sel_dinf_in_fsh(1,nfsh,1,nyrs)

  vector seldec_pen_fsh(1,nfsh) ;
  vector nnodes_fsh(1,nfsh) ;
  int seldecage ;
  !! seldecage = int(nages/2);
  ivector nselages_in_fsh(1,nfsh)

  ivector n_sel_ch_fsh(1,nfsh);
  ivector n_sel_ch_ind(1,nind);
  imatrix yrs_sel_ch_tmp(1,nind,1,endyr-styr+1);
  imatrix yrs_sel_ch_tmp_ind(1,nind,1,endyr-styr+1);
  !! yrs_sel_ch_tmp_ind.initialize();

  ivector   ind_sel_opt(1,nind)
  ivector phase_sel_ind(1,nind)

  vector   curv_pen_ind(1,nind)

  matrix   logsel_slp_in_ind(1,nind,1,nyrs)
  matrix   sel_inf_in_ind(1,nind,1,nyrs)
  matrix   sel_dslp_in_ind(1,nind,1,nyrs)
  matrix   logsel_dslp_in_ind(1,nind,1,nyrs)
  matrix   sel_dinf_in_ind(1,nind,1,nyrs)
  matrix   sel_slp_in_ind(1,nind,1,nyrs)

  vector   logsel_slp_in_indv(1,nind)
  vector   sel_inf_in_indv(1,nind)
  vector   logsel_dslp_in_indv(1,nind)
  vector   sel_dinf_in_indv(1,nind)


  vector seldec_pen_ind(1,nind) ;
  matrix sel_change_in_ind(1,nind,styr,endyr);
  ivector nselages_in_ind(1,nind)
  matrix sel_change_in_fsh(1,nfsh,styr,endyr);
  imatrix yrs_sel_ch_fsh(1,nfsh,1,endyr-styr);
  matrix sel_sigma_fsh(1,nfsh,1,endyr-styr);
  imatrix yrs_sel_ch_ind(1,nind,1,endyr-styr);
  matrix sel_sigma_ind(1,nind,1,endyr-styr);
  !! yrs_sel_ch_fsh.initialize();
  !! yrs_sel_ch_ind.initialize();
  !! sel_sigma_fsh.initialize();
  !! sel_sigma_ind.initialize();

  // Phase of estimation
  ivector phase_selcoff_fsh(1,nfsh)
  ivector phase_logist_fsh(1,nfsh)
  ivector phase_dlogist_fsh(1,nfsh)
  ivector phase_sel_spl_fsh(1,nfsh)

  ivector phase_selcoff_ind(1,nind)
  ivector phase_logist_ind(1,nind)
  ivector phase_dlogist_ind(1,nind)
  vector  sel_fsh_tmp(1,nages); 
  vector  sel_ind_tmp(1,nages); 
  3darray log_selcoffs_fsh_in(1,nfsh,1,nyrs,1,nages)
  3darray log_selcoffs_ind_in(1,nind,1,nyrs,1,nages)
  3darray  log_sel_spl_fsh_in(1,nfsh,1,nyrs,1,nages) // use nages for input to start
  // 3darray log_selcoffs_ind_in(1,nind,1,nyrs,1,nages)

 LOCAL_CALCS
  logsel_slp_in_fshv.initialize();
  sel_inf_in_fshv.initialize();
  logsel_dslp_in_fshv.initialize();
  sel_inf_in_fshv.initialize();
  sel_dinf_in_fshv.initialize();

  sel_inf_in_indv.initialize();
  logsel_dslp_in_indv.initialize();
  sel_inf_in_indv.initialize();
  sel_dinf_in_indv.initialize();

  phase_selcoff_ind.initialize();
  phase_logist_ind.initialize();
  phase_dlogist_ind.initialize();
  sel_fsh_tmp.initialize() ;
  sel_ind_tmp.initialize() ;
  log_selcoffs_fsh_in.initialize();
  log_selcoffs_ind_in.initialize();

  // nselages_in_fsh.initialize()   ;  
  // nselages_in_ind.initialize()   ;  
  nselages_in_fsh = nages-1;
  nselages_in_ind = nages-1;
  sel_change_in_ind.initialize()   ;  
  sel_slp_in_fsh.initialize()   ;  // ji
  sel_slp_in_ind.initialize()   ;  // ji
  sel_inf_in_fsh.initialize()   ;  // ji
  sel_inf_in_ind.initialize()   ;  // ji
  logsel_slp_in_fsh.initialize();  // ji
  logsel_slp_in_fshv.initialize();  // ji
  logsel_dslp_in_fsh.initialize(); // ji
  logsel_slp_in_ind.initialize();  // ji
  logsel_slp_in_indv.initialize();  // ji
  logsel_dslp_in_ind.initialize(); // ji
  sel_change_in_fsh.initialize()   ;  
  for (k=1;k<=nfsh;k++)
  {
    *(ad_comm::global_datafile) >> fsh_sel_opt(k)  ;  
    log_input(fsh_sel_opt(k));
    switch (fsh_sel_opt(k))
    {
      case 1 : // Selectivity coefficients 
      {
        *(ad_comm::global_datafile) >> nselages_in_fsh(k)   ;  
        *(ad_comm::global_datafile) >> phase_sel_fsh(k);  
        *(ad_comm::global_datafile) >> curv_pen_fsh(k) ;
        *(ad_comm::global_datafile) >> seldec_pen_fsh(k) ;
        seldec_pen_fsh(k) *= seldec_pen_fsh(k) ;  // square the input penalty
        *(ad_comm::global_datafile) >>  n_sel_ch_fsh(k) ;  
        n_sel_ch_fsh(k) +=1;
        yrs_sel_ch_fsh(k,1) = styr; // first year always estimated
        for (int i=2;i<=n_sel_ch_fsh(k);i++)
          *(ad_comm::global_datafile) >>  yrs_sel_ch_fsh(k,i) ;  
        for (int i=2;i<=n_sel_ch_fsh(k);i++)
          *(ad_comm::global_datafile) >>  sel_sigma_fsh(k,i) ;  
        log_input(nselages_in_fsh(k)) ;  
        log_input(phase_sel_fsh(k)) ;  
        log_input(curv_pen_fsh(k)) ;  
        log_input(seldec_pen_fsh(k)) ;  
        log_input(n_sel_ch_fsh(k)) ;  
        log_input(yrs_sel_ch_fsh(k)) ;  
        log_input(sel_sigma_fsh(k)) ;  
        // for (int i=styr;i<=endyr;i++) *(ad_comm::global_datafile) >> sel_change_in_fsh(k,i) ;
        sel_change_in_fsh(k,styr)=1.; 
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        // This to read in pre-specified selectivity values...
        sel_fsh_tmp.initialize();
        log_selcoffs_fsh_in.initialize();
        for (int j=1;j<=nages;j++) 
          *(ad_comm::global_datafile) >> sel_fsh_tmp(j);  
        log_selcoffs_fsh_in(k,1)(1,nselages_in_fsh(k)) = log((sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7)/mean(sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7) );
        for (int jj=2;jj<=n_sel_ch_fsh(k);jj++) 
        {
          // Set the selectivity for the oldest group
          for (int j=nselages_in_fsh(k)+1;j<=nages;j++) 
          {
            sel_fsh_tmp(j)  = sel_fsh_tmp(nselages_in_fsh(k));  
          }
          // Set tmp to actual initial vectors...
          log_selcoffs_fsh_in(k,jj)(1,nselages_in_fsh(k)) = log((sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7)/mean(sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7) );
          write_input_log<<"Sel_in_fsh "<< mfexp(log_selcoffs_fsh_in(k,jj))<<endl;
        }
        // exit(1);
        phase_selcoff_fsh(k) = phase_sel_fsh(k);
        phase_logist_fsh(k)  = -1;
        phase_dlogist_fsh(k) = -1;
        phase_sel_spl_fsh(k) = -1;
      }
        break;
      case 2 : // Single logistic
      {
        *(ad_comm::global_datafile) >> phase_sel_fsh(k);  
        *(ad_comm::global_datafile) >>  n_sel_ch_fsh(k) ;  
        n_sel_ch_fsh(k) +=1;
        yrs_sel_ch_fsh(k,1) = styr;
        for (int i=2;i<=n_sel_ch_fsh(k);i++)
          *(ad_comm::global_datafile) >>  yrs_sel_ch_fsh(k,i) ;  
        for (int i=2;i<=n_sel_ch_fsh(k);i++)
          *(ad_comm::global_datafile) >>  sel_sigma_fsh(k,i) ;  
        // This to read in pre-specified selectivity values...
        *(ad_comm::global_datafile) >> sel_slp_in_fsh(k,1) ;
        *(ad_comm::global_datafile) >> sel_inf_in_fsh(k,1) ;
        logsel_slp_in_fsh(k,1)   = log(sel_slp_in_fsh(k,1)) ;
        for (int jj=2;jj<=n_sel_ch_fsh(k);jj++) 
        {
          sel_inf_in_fsh(k,jj)    =     sel_inf_in_fsh(k,1) ;
          logsel_slp_in_fsh(k,jj) = log(sel_slp_in_fsh(k,1)) ;
        }
        log_input(phase_sel_fsh(k));
        log_input(n_sel_ch_fsh(k));
        log_input(sel_slp_in_fsh(k)(1,n_sel_ch_fsh(k)));
        log_input(sel_inf_in_fsh(k)(1,n_sel_ch_fsh(k)));
        log_input(logsel_slp_in_fsh(k)(1,n_sel_ch_fsh(k)));
        log_input(yrs_sel_ch_fsh(k)(1,n_sel_ch_fsh(k)));

        phase_selcoff_fsh(k) = -1;
        phase_logist_fsh(k) = phase_sel_fsh(k);
        phase_dlogist_fsh(k) = -1;
        phase_sel_spl_fsh(k) = -1;

        logsel_slp_in_fshv(k) = logsel_slp_in_fsh(k,1);
           sel_inf_in_fshv(k) =    sel_inf_in_fsh(k,1);
        break;
      }
      case 3 : // Double logistic 
      {
        write_input_log << "Double logistic abandoned..."<<endl;exit(1);
        break;
      }
      case 4 : // Splines         
      {
      }
      break;
      write_input_log << fshname(k)<<" fish sel opt "<<endl<<fsh_sel_opt(k)<<" "<<endl<<"Sel_change"<<endl<<sel_change_in_fsh(k)<<endl;
    }
  }
  // Indices here..............
  yrs_sel_ch_ind.initialize() ;  
  sel_sigma_ind.initialize();
  for(k=1;k<=nind;k++)
  {
    *(ad_comm::global_datafile) >> ind_sel_opt(k)  ;  
    write_input_log << endl<<"Survey "<<indname(k)<<endl;
    log_input(ind_sel_opt(k));
    switch (ind_sel_opt(k))
    {
      case 1 : // Selectivity coefficients  indices
      {
        *(ad_comm::global_datafile) >> nselages_in_ind(k)   ;  
        *(ad_comm::global_datafile) >> phase_sel_ind(k);  
        *(ad_comm::global_datafile) >> curv_pen_ind(k) ;
        *(ad_comm::global_datafile) >> seldec_pen_ind(k) ;
        seldec_pen_ind(k) *= seldec_pen_ind(k);
        *(ad_comm::global_datafile) >>  n_sel_ch_ind(k) ;  
        n_sel_ch_ind(k)+=1;
        yrs_sel_ch_ind(k,1) = styr;
        yrs_sel_ch_tmp_ind(k,1) = styr;
        for (int i=2;i<=n_sel_ch_ind(k);i++)
          *(ad_comm::global_datafile) >>  yrs_sel_ch_ind(k,i) ;  
        for (int i=2;i<=n_sel_ch_ind(k);i++)
          *(ad_comm::global_datafile) >>  sel_sigma_ind(k,i) ;  
        sel_change_in_ind(k,styr)=1.; 
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        log_input(indname(k));
        log_input(nselages_in_ind(k));
        log_input(phase_sel_ind(k));
        log_input(seldec_pen_ind(k));
        log_input(n_sel_ch_ind(k));
        log_input(sel_change_in_ind(k));
        log_input(n_sel_ch_ind(k));
        // log_input(yrs_sel_ch_ind(k)(1,n_sel_ch_ind(k)));
        log_input(yrs_sel_ch_ind(k));
        // This to read in pre-specified selectivity values...
        for (j=1;j<=nages;j++) 
          *(ad_comm::global_datafile) >> sel_ind_tmp(j);  
        log_input(sel_ind_tmp);
        log_selcoffs_ind_in(k,1)(1,nselages_in_ind(k)) = log((sel_ind_tmp(1,nselages_in_ind(k))+1e-7)/mean(sel_ind_tmp(1,nselages_in_ind(k))+1e-7) );
        // set all change selectivity to initial values
        for (int jj=2;jj<=n_sel_ch_ind(k);jj++) 
        {
          for (int j=nselages_in_ind(k)+1;j<=nages;j++) // This might be going out of nages=nselages
          {
            sel_ind_tmp(j)  = sel_ind_tmp(nselages_in_ind(k));  
          }
          // Set tmp to actual initial vectors...
          log_selcoffs_ind_in(k,jj)(1,nselages_in_ind(k)) = log((sel_ind_tmp(1,nselages_in_ind(k))+1e-7)/mean(sel_ind_tmp(1,nselages_in_ind(k))+1e-7) );
          write_input_log<<"Sel_in_ind "<< mfexp(log_selcoffs_ind_in(k,jj))<<endl;
        }
        phase_selcoff_ind(k) = phase_sel_ind(k);
        phase_logist_ind(k)  = -2;
        phase_dlogist_ind(k) = -1;
      }
      break;
      case 2 : // Single logistic
      {
        *(ad_comm::global_datafile) >> phase_sel_ind(k);  
        *(ad_comm::global_datafile) >>  n_sel_ch_ind(k) ;  
        n_sel_ch_ind(k) +=1;
        yrs_sel_ch_ind(k,1) = styr; // first year always estimated
        yrs_sel_ch_tmp_ind(k,1) = styr;
        for (int i=2;i<=n_sel_ch_ind(k);i++)
          *(ad_comm::global_datafile) >>  yrs_sel_ch_ind(k,i) ;  
        for (int i=2;i<=n_sel_ch_ind(k);i++)
          *(ad_comm::global_datafile) >>  sel_sigma_ind(k,i) ;  
        sel_change_in_ind(k,styr)=1.; 

        log_input(indname(k));
        log_input(nselages_in_ind(k));
        log_input(phase_sel_ind(k));
        log_input(sel_change_in_ind(k));
        log_input(n_sel_ch_ind(k));
        log_input(yrs_sel_ch_ind(k)(1,n_sel_ch_ind(k)));
        // This to read in pre-specified selectivity values...
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) { if(sel_change_in_ind(k,i)>0) { j++; yrs_sel_ch_tmp_ind(k,j) = i; } }
        // This to read in pre-specified selectivity values...
        *(ad_comm::global_datafile) >> sel_slp_in_ind(k,1) ;
        *(ad_comm::global_datafile) >> sel_inf_in_ind(k,1) ;
        logsel_slp_in_ind(k,1) =   log(sel_slp_in_ind(k,1)) ;
        for (int jj=2;jj<=n_sel_ch_ind(k);jj++) 
        {
          sel_inf_in_ind(k,jj)    =     sel_inf_in_ind(k,1) ;
          logsel_slp_in_ind(k,jj) = log(sel_slp_in_ind(k,1)) ;
        }
        log_input(sel_slp_in_ind(k,1));
        log_input(sel_inf_in_ind(k,1));
        log_input(logsel_slp_in_ind(k,1));

        phase_selcoff_ind(k) = -1;
        phase_logist_ind(k) = phase_sel_ind(k);
        phase_dlogist_ind(k)  = -1;

        logsel_slp_in_indv(k) = logsel_slp_in_ind(k,1);
           sel_inf_in_indv(k) =    sel_inf_in_ind(k,1);
        log_input(logsel_slp_in_indv(k));
      }
      break;
      case 3 : // Double logistic 
      {
        write_input_log << "Double logistic abandoned..."<<endl;exit(1);
      }
        break;
      case 4 : // spline for indices
      {
      }
      break;
    }
    write_input_log << indname(k)<<" ind sel opt "<<ind_sel_opt(k)<<" "<<sel_change_in_ind(k)<<endl;
  }
  write_input_log<<"Phase indices Sel_Coffs: "<<phase_selcoff_ind<<endl; 
 END_CALCS
  init_matrix wt_pop(1,nstk,1,nages)
  !! log_input(wt_pop);
  init_matrix maturity(1,nstk,1,nages)
  !! log_input(maturity);
  // !! if (max(maturity)>.9) maturity /=2.;
  matrix wt_mature(1,nstk,1,nages);
  !! for (s=1;s<=nstk;s++)
  !!  wt_mature(s) = elem_prod(wt_pop(s),maturity(s)) ;
  init_number test;
  !! write_input_log<<" Test: "<<test<<endl;
 !! if (test!=123456789) {cerr<<"Control file not read in correctly... "<<endl;exit(1);}


  ivector nopt_fsh(1,2) // number of options...
  !! nopt_fsh.initialize();
  !! for (k=1;k<=nfsh;k++) if(fsh_sel_opt(k)==1) nopt_fsh(1)++;else nopt_fsh(2)++;

  // Fishery selectivity description:
  // type 1
  
  // Number of ages

  !! write_input_log << "# Fshry Selages: " << nselages_in_fsh  <<endl;
  !! write_input_log << "# Srvy  Selages: " << nselages_in_ind <<endl;



  !! write_input_log << "# Phase for age-spec fishery "<<phase_selcoff_fsh<<endl;
  !! write_input_log << "# Phase for logistic fishery "<<phase_logist_fsh<<endl;
  !! write_input_log << "# Phase for dble logistic fishery "<<phase_dlogist_fsh<<endl;

  !! write_input_log << "# Phase for age-spec indices  "<<phase_selcoff_ind<<endl;
  !! write_input_log << "# Phase for logistic indices  "<<phase_logist_ind<<endl;
  !! write_input_log << "# Phase for dble logistic ind "<<phase_dlogist_ind<<endl;

  !! for (k=1;k<=nfsh;k++) if (phase_selcoff_fsh(k)>0) curv_pen_fsh(k) = 1./ (square(curv_pen_fsh(k))*2.);
  !! write_input_log<<"# Curv_pen_fsh: "<<endl<<curv_pen_fsh<<endl;
  !! for (k=1;k<=nind;k++) if (phase_selcoff_ind(k)>0) curv_pen_ind(k) = 1./ (square(curv_pen_ind(k))*2.);
  !! write_input_log<<"# Curv_pen_ind: "<<endl<<curv_pen_fsh<<endl;

  int  phase_fmort;
  int  phase_proj;
  ivector   nselages_fsh(1,nfsh);
  matrix xnodes_fsh(1,nfsh,1,nnodes_fsh)
  matrix xages_fsh(1,nfsh,1,nages)

  ivector   nselages_ind(1,nind);
  //Resetting data here for retrospectives////////////////////////////////////////////
 LOCAL_CALCS
  for (int k=1;k<=nfsh;k++) 
  {
    // if ((endyr-retro)<=yrs_sel_ch_fsh(k,n_sel_ch_fsh(k))) n_sel_ch_fsh(k)-=retro ;  //Ojo not always true
    for (int i=1;i<=retro;i++)
    {
      if (n_sel_ch_fsh(k) >0)
        if (yrs_sel_ch_fsh(k,n_sel_ch_fsh(k))>=(endyr-retro))
          n_sel_ch_fsh(k) -= 1;
    }
    for (int i=1;i<=retro;i++) 
    {
      if (nyrs_fsh_age(k) >0)
        if (max(yrs_fsh_age_in(k)(1,nyrs_fsh_age(k)))>=(endyr-retro)) 
          nyrs_fsh_age(k) -= 1;
      if (nyrs_fsh_length(k) >0)
        if (max(yrs_fsh_length_in(k)(1,nyrs_fsh_length(k)))>=(endyr-retro)) 
          nyrs_fsh_length(k) -= 1;
    }
  }
  // now for indices
  for (int k=1;k<=nind;k++) 
  {
    // if ((endyr-retro)<=yrs_sel_ch_ind(k,n_sel_ch_ind(k))) n_sel_ch_ind(k)-=retro ;  //Ojo not always true
    for (int i=1;i<=retro;i++)
    {
      if (n_sel_ch_ind(k) >0)
        if (yrs_sel_ch_ind(k,n_sel_ch_ind(k))>=(endyr-retro)) 
          n_sel_ch_ind(k) -= 1;
    }
    for (int i=1;i<=retro;i++) 
    {
      // index values
      if (nyrs_ind(k) >0)
        if (max(yrs_ind_in(k)(1,nyrs_ind(k)))>=(endyr-retro)) 
          nyrs_ind(k) -= 1;
      // Ages (since they can be different than actual index years)
      if (nyrs_ind_age(k) >0)
        if (max(yrs_ind_age_in(k)(1,nyrs_ind_age(k)))>=(endyr-retro)) 
          nyrs_ind_age(k) -= 1;
      // Lengths
      if (nyrs_ind_length(k) >0)
        if (max(yrs_ind_length_in(k)(1,nyrs_ind_length(k)))>=(endyr-retro)) 
          nyrs_ind_length(k) -= 1;
    }
  }
  for (s=1;s<=nstk;s++)
  {
    endyr_rec_est(s,nreg(s)) = endyr_rec_est(s,nreg(s)) - retro;
  }
  endyr         = endyr - retro;
  styr_fut      = endyr+1;
  endyr_fut     = endyr + nproj_yrs; 
  endyr_sp      = endyr   - rec_age - 1;// endyr year of (main) spawning biomass
  log_input(styr_fut);
  log_input(endyr_fut);
  log_input(nyrs_fsh_age);
 END_CALCS
 // now use redimensioned data for retro
  matrix catch_bio(1,nfsh,styr,endyr)         //Catch biomass 
  matrix catch_bio_sd(1,nfsh,styr,endyr)      //Catch biomass standard errors 
  matrix catch_bio_lsd(1,nfsh,styr,endyr)     //Catch biomass standard errors (for lognormal)
  matrix catch_bio_lva(1,nfsh,styr,endyr)     //Catch biomass variance (for lognormal)
  matrix catch_bioT(styr,endyr,1,nfsh)
  vector catch_lastyr(1,nfsh);
  imatrix yrs_fsh_age(1,nfsh,1,nyrs_fsh_age)
  imatrix yrs_fsh_length(1,nfsh,1,nyrs_fsh_length)
  matrix  n_sample_fsh_age(1,nfsh,1,nyrs_fsh_age)    //Years of index index value (annual)
  matrix n_sample_fsh_length(1,nfsh,1,nyrs_fsh_length)    //Years of index index value (annual)
  3darray oac_fsh(1,nfsh,1,nyrs_fsh_age,1,nages)
  3darray olc_fsh(1,nfsh,1,nyrs_fsh_length,1,nlength)

  imatrix yrs_ind(1,nind,1,nyrs_ind)         //Years of index value (annual)
  matrix obs_ind(1,nind,1,nyrs_ind)          //values of index value (annual)
  matrix obs_se_ind(1,nind,1,nyrs_ind)       //values of indices serrs

  imatrix yrs_ind_age(1,nind,1,nyrs_ind_age)  //Years of index value (annual)
  imatrix yrs_ind_length(1,nind,1,nyrs_ind_length)
  matrix n_sample_ind_age(1,nind,1,nyrs_ind_age)         //Years of index value (annual)
  matrix n_sample_ind_length(1,nind,1,nyrs_ind_length)    //Years of index index value (annual)
  3darray oac_ind(1,nind,1,nyrs_ind_age,1,nages)  //values of Index proportions at age
  3darray olc_ind(1,nind,1,nyrs_ind_length,1,nlength)

  matrix     obs_lse_ind(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
  matrix     obs_lva_ind(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
 LOCAL_CALCS
  for (s=1;s<=nstk;s++)
  {
    yy_shift_st(s,1) = styr_rec;
    yy_shift_end(s,nreg(s)) = endyr;
    for (r=2;r<=nreg(s);r++)
    {
      yy_shift_st(s,r) = reg_shift(s,r-1);
      yy_shift_end(s,r-1) = reg_shift(s,r-1)-1;
    }
  }
  for (int k=1;k<=nfsh;k++)
  {
    catch_bio(k) = catch_bio_in(k)(styr,endyr);
    catch_bio_sd(k) = catch_bio_sd_in(k)(styr,endyr);
    if (nyrs_fsh_age(k)>0)
    {
      yrs_fsh_age(k) = yrs_fsh_age_in(k)(1,nyrs_fsh_age(k));
      n_sample_fsh_age(k) = n_sample_fsh_age_in(k)(1,nyrs_fsh_age(k));
    }
    if (nyrs_fsh_length(k)>0)
    {
      yrs_fsh_length(k) = yrs_fsh_length_in(k)(1,nyrs_fsh_length(k));
      n_sample_fsh_length(k) = n_sample_fsh_length_in(k)(1,nyrs_fsh_length(k));
    }
    for (int i=1;i<=nyrs_fsh_age(k);i++)
      oac_fsh(k,i) = oac_fsh_in(k,i) ;
    for (int i=1;i<=nyrs_fsh_length(k);i++)
      olc_fsh(k,i) = olc_fsh_in(k,i) ;
  }
  catch_bio_lsd = sqrt(log(square(catch_bio_sd) + 1.));
  catch_bio_lva = log(square(catch_bio_sd) + 1.);
  catch_bioT    = trans(catch_bio);
  catch_lastyr  = catch_bioT(endyr);
  for (int k=1;k<=nind;k++)
  {
    yrs_ind(k)  = yrs_ind_in(k)(1,nyrs_ind(k));
    obs_ind(k)  = obs_ind_in(k)(1,nyrs_ind(k));
    obs_se_ind(k)  = obs_se_ind_in(k)(1,nyrs_ind(k));

    if (nyrs_ind_age(k)>0)
    {
      yrs_ind_age(k) = yrs_ind_age_in(k)(1,nyrs_ind_age(k));
      n_sample_ind_age(k) = n_sample_ind_age_in(k)(1,nyrs_ind_age(k));
    }
    if (nyrs_ind_length(k)>0)
    {
      yrs_ind_length(k) = yrs_ind_length_in(k)(1,nyrs_ind_length(k));
      n_sample_ind_length(k) = n_sample_ind_length_in(k)(1,nyrs_ind_length(k));
  }
    for (int i=1;i<=nyrs_ind_age(k);i++)
      oac_ind(k,i) = oac_ind_in(k,i) ;
    for (int i=1;i<=nyrs_ind_length(k);i++)
      olc_ind(k,i) = olc_ind_in(k,i) ;
  }
  log_input(nyrs_fsh_age);
  log_input(yrs_fsh_age);
  log_input(n_sample_fsh_age);
  log_input(oac_fsh);
  log_input(olc_fsh);
  log_input(wt_fsh);

  log_input(nyrs_ind_age);
  log_input(yrs_ind_age);
  log_input(n_sample_ind_age);
  log_input(oac_ind);
  log_input(olc_ind);
  obs_lse_ind = elem_div(obs_se_ind,obs_ind);
  obs_lse_ind = sqrt(log(square(obs_lse_ind) + 1.));
  log_input(obs_lse_ind);
  obs_lva_ind = square(obs_lse_ind);
 END_CALCS

  ////////////////////////////////////////////////////////////////////////////////////
 LOCAL_CALCS
  for (k=1; k<=nfsh;k++)
  {
    // xages_fsh increments from 0-1 by number of ages, say
    xages_fsh.initialize();
    log_input(xages_fsh);
    xages_fsh(k).fill_seqadd(0.,1.0/(nages-1));
    log_input(xages_fsh);
    //  xnodes increments from 0-1 by number of nodes
    xnodes_fsh.initialize();
    xnodes_fsh(k).fill_seqadd(0.,1.0/(nnodes_fsh(k)-1));
    log_input(xnodes_fsh);
    // xages_fsh(k).fill_seqadd(0,1.0/(nselages_in_fsh(k)-1)); //prefer to use nselages but need 3d version to work
  }
  write_input_log<<"Yrs fsh_sel change: "<<yrs_sel_ch_fsh<<endl;
  // for (k=1; k<=nind;k++) yrs_sel_ch_ind(k) = yrs_sel_ch_tmp_ind(k)(1,n_sel_ch_ind(k));
  write_input_log<<"Yrs ind_sel change: "<<yrs_sel_ch_ind<<endl;
    log_sigmarprior = log(sigmarprior);
    log_input(steepnessprior);
    log_input(sigmarprior);
    for (s=1;s<=nstk;s++)
    {
      nrecs_est_shift(cum_regs(s)+nreg(s)) = endyr_rec_est(s,nreg(s))-styr_rec_est(s,nreg(s))+1;
    }
    write_input_log<<"#  SSB estimated in styr endyr: " <<styr_sp    <<" "<<endyr_sp      <<" "<<endl;
    write_input_log<<"#  Rec estimated in styr endyr: " <<styr_rec    <<" "<<endyr        <<" "<<endl;
    write_input_log<<"#  SR Curve fit  in styr endyr: " <<styr_rec_est<<endl<<" "<<endyr_rec_est<<" "<<endl;
    write_input_log<<"#             Model styr endyr: " <<styr        <<" "<<endyr        <<" "<<endl;
    log_qprior = log(qprior);
    log_input(qprior);
    log_q_power_prior = log(q_power_prior);
    write_input_log<<"# q_power_prior " <<endl<<q_power_prior<<" "<<endl;
    write_input_log<<"# cv_catchbiomass " <<endl<<cv_catchbiomass<<" "<<endl;
    write_input_log<<"# CatchbiomassPen " <<endl<<catchbiomass_pen<<" "<<endl;
    write_input_log<<"# Number of projection years " <<endl<<nproj_yrs<<" "<<endl;// cin>>junk;

 END_CALCS

  vector offset_ind(1,nind)
  vector offset_fsh(1,nfsh)
  vector offset_lfsh(1,nfsh)
  vector offset_lind(1,nind)

  int do_fmort;
  !! do_fmort=0;
  int Popes;
 LOCAL_CALCS
  Popes=0; // option to do Pope's approximation (not presently flagged outside of code)
  if (Popes) 
    phase_fmort = -2;
  else
    phase_fmort = 1;

  phase_proj  =  5;

  Steepness_UB = .9999; // upper bound of steepness
  offset_ind.initialize();
  offset_fsh.initialize();
  offset_lfsh.initialize();
  offset_lind.initialize();
  double sumtmp;
  for (k=1;k<=nfsh;k++)
    for (i=1;i<=nyrs_fsh_age(k);i++)
    {
      oac_fsh(k,i) /= sum(oac_fsh(k,i)); // Normalize to sum to one
      offset_fsh(k) -= n_sample_fsh_age(k,i)*(oac_fsh(k,i) + 0.001) * log(oac_fsh(k,i) + 0.001 ) ;
    }
  for (k=1;k<=nfsh;k++)
    for (i=1;i<=nyrs_fsh_length(k);i++)
    {
      olc_fsh(k,i) /= sum(olc_fsh(k,i)); // Normalize to sum to one
      offset_lfsh(k) -= n_sample_fsh_length(k,i)*(olc_fsh(k,i) + 0.001) * log(olc_fsh(k,i) + 0.001 ) ;
    }

  for (k=1;k<=nind;k++)
  {
    for (i=1;i<=nyrs_ind_age(k);i++)
    {
      oac_ind(k,i) /= sum(oac_ind(k,i)); // Normalize to sum to one
      offset_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(oac_ind(k,i) + 0.001 ) ;
    }
    for (i=1;i<=nyrs_ind_length(k);i++)
    {
      olc_ind(k,i) /= sum(olc_ind(k,i)); // Normalize to sum to one
      offset_lind(k) -= n_sample_ind_length(k,i)*(olc_ind(k,i) + 0.001) * log(olc_ind(k,i) + 0.001 ) ;
    }
  }
  log_input(offset_fsh); 
  log_input(offset_ind); 

  if (ad_comm::argc > 1) // Command line argument to profile Fishing mortality rates...
  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-uFmort"))>-1)
      do_fmort=1;
  }

  // Compute an initial Rzero value based on exploitation 
  for (s=1;s<=nstk;s++)
  {
    for (r=1;r<=nreg(s);r++)
    {
      double btmp=0.;
      double ctmp=0.;
      dvector ntmp(1,nages);
      ntmp(1) = 1.;
      for (int a=2;a<=nages;a++)
        ntmp(a) = ntmp(a-1)*exp(-natmortprior(mort_map(s,r))-.05);
      btmp = wt_pop(s) * ntmp;
      write_input_log << "Mean Catch "<< "stock "<< s << "regime "<< r <<endl;
      int yy_shift_st_tmp;
      int yy_shift_end_tmp;
      if (r > 1)
        yy_shift_st_tmp = yy_shift_st(s,r);
      else
        yy_shift_st_tmp = styr;
      yy_shift_end_tmp = yy_shift_end(s,r);
      int jj=0;
      for (k=1;k<=nfsh;k++)
      {
        if (sel_map(1,k) == s)
        {
          ctmp += sum(catch_bio(k)(yy_shift_st_tmp,yy_shift_end_tmp));
          jj++;
        }
      }
      ctmp /= (jj * (yy_shift_end_tmp - yy_shift_st_tmp + 1));
      write_input_log << ctmp <<endl;
     if( phase_mean_rec_tmp != -1 | s != istk_tmp | r != ireg_tmp)
       R_guess(cum_regs(s)+r) = log((ctmp/.02 )/btmp) ;
    }
  }
  write_input_log << "R_guess "<<endl;
  write_input_log << R_guess <<endl;
 END_CALCS
 // vector len_bins(1,nlength)
 // !! len_bins.fill_seqadd(stlength,binlength);

PARAMETER_SECTION
 // Biological Parameters
  init_bounded_number_vector Mest(1,nmort,.02,4.8,phase_M)
  init_bounded_vector_vector Mage_offset(1,nmort,1,npars_Mage,-3,3,phase_Mage)
  matrix Mage(1,nmort,1,nages)
  init_bounded_vector_vector M_rw(1,nstk,1,npars_rw_M,-10,10,phase_rw_M)
  matrix natmort(1,nstk,styr,endyr)
  3darray natage(1,nstk,styr,endyr+1,1,nages)
  3darray N_NoFsh(1,nstk,styr,endyr_fut,1,nages);
  // matrix Sp_Biom(1,nstk,styr_sp,endyr)
  matrix pred_rec(1,nstk,styr_rec,endyr)
  matrix mod_rec(1,nstk,styr_rec,endyr) // As estimated by model
  3darray  M(1,nstk,styr,endyr,1,nages)
  3darray  Z(1,nstk,styr,endyr,1,nages)
  3darray  S(1,nstk,styr,endyr,1,nages)


 //-----GROWTH PARAMETERS--------------------------------------------------
  init_number_vector log_Linf(1,ngrowth,phase_Linf);
  init_number_vector log_k(1,ngrowth,phase_k);
  init_number_vector log_Lo(1,ngrowth,phase_Lo);
  init_number_vector log_sdage(1,ngrowth,phase_sdage);
  matrix             wt_age_vb(1,ngrowth,1,nages)
  matrix             maturity_vb(1,ngrowth,1,nages)
//---------------------------------------------------------------------------


 // Stock rectuitment params
  init_number_vector mean_log_rec(1,nregs,phase_mean_rec); 
  init_bounded_number_vector steepness(1,nrec,0.21,Steepness_UB,phase_srec)
  init_number_vector log_Rzero(1,nregs,phase_Rzero)  
  // OjO
  // init_bounded_vector initage_dev(2,nages,-15,15,4)
  init_bounded_matrix rec_dev(1,nstk,styr_rec,endyr,-15,15,2)
  // init_matrix rec_dev(1,nstk,styr_rec,endyr,2)
  init_number_vector log_sigmar(1,nrec,phase_sigmar);
  vector m_sigmarsq(1,nregs)
  vector m_sigmar(1,nregs)
  vector sigmarsq(1,nregs)
  vector sigmar(1,nrec)
  vector alpha(1,nregs)
  vector beta(1,nregs)
  vector Bzero(1,nregs)
  vector Rzero(1,nregs)
  vector phizero(1,nregs)
  vector avg_rec_dev(1,nregs)

 // Fishing mortality parameters
  // init_vector         log_avg_fmort(1,nfsh,phase_fmort)
  // init_bounded_matrix fmort_dev(1,nfsh,styr,endyr,-15,15.,phase_fmort)
  init_bounded_matrix fmort(1,nfsh,styr,endyr,0.00,5.,phase_fmort)
  matrix Fmort(1,nstk,styr,endyr);  // Annual total Fmort
  number hrate
  number catch_tmp
  number Fnew 

  !! for (k=1;k<=nfsh;k++) nselages_fsh(k)=nselages_in_fsh(k); // Sets all elements of a vector to one scalar value...
  !! for (k=1;k<=nind;k++) nselages_ind(k)=nselages_in_ind(k); // Sets all elements of a vector to one scalar value...

 //  init_3darray log_selcoffs_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_selcoff_fsh)
  init_matrix_vector log_selcoffs_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_selcoff_fsh) // 3rd dimension out...
  // option to estimate smoother for selectivity penalty
  // init_number_vector logSdsmu_fsh(1,nfsh,1,phase_selcoff_fsh) 
  !! if (fsh_sel_opt(1)==4) nnodes_tmp=nnodes_fsh(1);  // NOTE THIS won't work in general
  //init_matrix_vector  log_sel_spl_fsh(1,nfsh,1,n_sel_ch_fsh,1,nnodes_tmp,phase_sel_spl_fsh)
  init_matrix_vector  log_sel_spl_fsh(1,nfsh,1,n_sel_ch_fsh,1,4,phase_sel_spl_fsh)

  !! log_input(nfsh);
  !! log_input(n_sel_ch_fsh);
  !! log_input(nselages_fsh);
  !! log_input(phase_selcoff_fsh);
  init_vector_vector logsel_slope_fsh(1,nfsh,1,n_sel_ch_fsh,phase_logist_fsh)
  matrix                sel_slope_fsh(1,nfsh,1,n_sel_ch_fsh)
  init_vector_vector     sel50_fsh(1,nfsh,1,n_sel_ch_fsh,phase_logist_fsh)
  init_vector_vector logsel_dslope_fsh(1,nfsh,1,n_sel_ch_fsh,phase_dlogist_fsh)
  matrix                sel_dslope_fsh(1,nfsh,1,n_sel_ch_fsh)
  !! int lb_d50=nages/2;
  init_bounded_vector_vector     seld50_fsh(1,nfsh,1,n_sel_ch_fsh,lb_d50,nages,phase_dlogist_fsh)

  // !!exit(1);
  3darray log_sel_fsh(1,nfsh,styr,endyr,1,nages)
  3darray sel_fsh(1,nfsh,styr,endyr,1,nages)
  matrix avgsel_fsh(1,nfsh,1,n_sel_ch_fsh);

  3darray Ftot(1,nstk,styr,endyr,1,nages)
  3darray F(1,nfsh,styr,endyr,1,nages)
  3darray eac_fsh(1,nfsh,1,nyrs_fsh_age,1,nages)
//-----------------------------------------------NEW--------
  3darray elc_fsh(1,nfsh,1,nyrs_fsh_length,1,nlength)
  3darray elc_ind(1,nind,1,nyrs_ind_length,1,nlength)
//----------------------------------------------------------
  matrix  pred_catch(1,nfsh,styr,endyr)
  3darray catage(1,nfsh,styr,endyr,1,nages)
  3darray catage_tot(1,nstk,styr,endyr,1,nages)
  matrix expl_biom(1,nfsh,styr,endyr)

 // Parameters for computing SPR rates 
  vector F50(1,nfsh)
  vector F40(1,nfsh)
  vector F35(1,nfsh)

 // Stuff for SPR and yield projections
  vector sigmar_fut(1,nstk)
  vector f_tmp(1,nfsh)
  number SB0
  number SBF50
  number SBF40
  number SBF35
  vector Fratio(1,nfsh)
  !! Fratio = 1;
  !! Fratio /= sum(Fratio);

  matrix Nspr(1,4,1,nages)
 
  3darray nage_future(1,nstk,styr_fut,endyr_fut,1,nages)

  init_matrix rec_dev_future(1,nstk,styr_fut,endyr_fut,phase_proj);
  matrix Sp_Biom_future(1,nstk,styr_fut-rec_age,endyr_fut);
  3darray F_future(1,nfsh,styr_fut,endyr_fut,1,nages);
  3darray Z_future(1,nstk,styr_fut,endyr_fut,1,nages);
  3darray S_future(1,nstk,styr_fut,endyr_fut,1,nages);
  3darray catage_future(1,nstk,styr_fut,endyr_fut,1,nages);
  vector avg_rec_dev_future(1,nstk)
  vector avg_F_future(1,5)

 // Survey Observation parameters
  init_number_vector log_q_ind(1,nind,phase_q) 
  init_number_vector log_q_power_ind(1,nind,phase_q_power) 
  init_vector_vector log_rw_q_ind(1,nind,1,npars_rw_q,phase_rw_q) 
  init_matrix_vector log_selcoffs_ind(1,nind,1,n_sel_ch_ind,1,nselages_ind,phase_selcoff_ind)

  // init_vector_vector logsel_slope_ind(1,nind,1,n_sel_ch_ind,phase_logist_ind) // Need to make positive or reparameterize
  init_vector_vector logsel_slope_ind(1,nind,1,n_sel_ch_ind,phase_logist_ind+1) // Need to make positive or reparameterize
  init_bounded_vector_vector        sel50_ind(1,nind,1,n_sel_ch_ind,1,20,phase_logist_ind)

  init_vector_vector  logsel_dslope_ind(1,nind,1,n_sel_ch_ind,phase_dlogist_ind) // Need to make positive or reparameterize
  init_bounded_vector_vector seld50_ind(1,nind,1,n_sel_ch_ind,lb_d50,nages,phase_dlogist_ind)

  matrix                sel_slope_ind(1,nind,1,n_sel_ch_ind)
  matrix                sel_dslope_ind(1,nind,1,n_sel_ch_ind)

  3darray log_sel_ind(1,nind,styr,endyr,1,nages)
  3darray sel_ind(1,nind,styr,endyr,1,nages)
  matrix avgsel_ind(1,nind,1,n_sel_ch_ind);

  matrix pred_ind(1,nind,1,nyrs_ind)
  3darray eac_ind(1,nind,1,nyrs_ind_age,1,nages)

 // Likelihood value names         
  number sigma
  matrix rec_like(1,nstk,1,4)
  vector catch_like(1,nfsh)
  vector age_like_fsh(1,nfsh)
//---------------------------------NEW
  vector length_like_fsh(1,nfsh)
  vector length_like_ind(1,nind)
//---------------------------------NEW

  vector age_like_ind(1,nind)
  matrix sel_like_fsh(1,nfsh,1,4)       
  matrix sel_like_ind(1,nind,1,4)       
  vector ind_like(1,nind)
  vector fpen(1,6)    
  matrix post_priors(1,npar,1,8)
  vector post_priors_indq(1,nind)
  objective_function_value obj_fun
  vector obj_comps(1,14)
  init_vector repl_F(1,nstk,5)

  sdreport_vector repl_yld(1,nstk)
  sdreport_vector repl_SSB(1,nstk)
  sdreport_vector B100(1,nstk)
  vector F50_est(1,nstk)
  vector F40_est(1,nstk)
  vector F35_est(1,nstk)
  matrix q_ind(1,nind,1,nyrs_ind)
  vector q_power_ind(1,nind)
  // sdreport_vector q_ind(1,nind)
  sdreport_matrix totbiom(1,nstk,styr,endyr+1)
  sdreport_matrix totbiom_NoFish(1,nstk,styr,endyr)
  sdreport_matrix Sp_Biom(1,nstk,styr_sp,endyr+1)
  sdreport_matrix Sp_Biom_NoFish(1,nstk,styr_sp,endyr_fut)
  sdreport_matrix Sp_Biom_NoFishRatio(1,nstk,styr,endyr)
  sdreport_vector ABCBiom(1,nstk);
  sdreport_matrix recruits(1,nstk,styr,endyr+1)
  // matrix recruits(1,nstk,styr,endyr+1)
  sdreport_vector depletion(1,nstk)
  sdreport_vector depletion_dyn(1,nstk)
  sdreport_vector MSY(1,nstk);
  sdreport_vector MSYL(1,nstk);
  sdreport_vector Fmsy(1,nstk);
  sdreport_vector lnFmsy(1,nstk);
  sdreport_vector Fcur_Fmsy(1,nstk);
  sdreport_vector Rmsy(1,nstk);
  sdreport_vector Bmsy(1,nstk);
  sdreport_vector Bcur_Bmsy(1,nstk);
  sdreport_vector pred_ind_nextyr(1,nind);
  sdreport_vector OFL(1,nstk);
  // NOTE TO DAVE: Need to have a phase switch for sdreport variables(
  3darray catch_future(1,nstk,1,4,styr_fut,endyr_fut); // Note, don't project for F=0 (it will bomb)
  //sdreport_matrix SSB_fut(1,5,styr_fut,endyr_fut) //Ojo
  3darray SSB_fut(1,nstk,1,5,styr_fut,endyr_fut) //Ojo
  sdreport_matrix SSB_fut_1(1,nstk,styr_fut,endyr_fut)
  sdreport_matrix SSB_fut_2(1,nstk,styr_fut,endyr_fut)
  sdreport_matrix SSB_fut_3(1,nstk,styr_fut,endyr_fut)
  sdreport_matrix SSB_fut_4(1,nstk,styr_fut,endyr_fut)
  sdreport_matrix SSB_fut_5(1,nstk,styr_fut,endyr_fut)
  !! write_input_log <<"logRzero "<<log_Rzero<<endl;
  !! write_input_log <<"logmeanrec "<<mean_log_rec<<endl;
  !! write_input_log<< "exp(log_sigmarprior "<<exp(log_sigmarprior)<<endl;




//-----GROWTH PARAMETERS--------------------------------------------------
 vector Linf(1,ngrowth);
 vector k_coeff(1,ngrowth);
 vector Lo(1,ngrowth);
 vector sdage(1,ngrowth);
 matrix mu_age(1,ngrowth,1,nages);
 matrix sigma_age(1,ngrowth,1,nages);
 matrix P1(1,nages,1,nlength);
 matrix P2(1,nages,1,nlength);
 matrix P3(1,nages,1,nlength);
 vector Ones_length(1,nlength);
 3darray P_age2len(1,ngrowth,1,nages,1,nlength);

//-----------------------------------------------------------------------
 // Initialize coefficients (if needed)
 LOCAL_CALCS
  for (k=1;k<=nfsh;k++) 
  {
    write_input_log<<"Fish sel phase: "<<phase_selcoff_fsh(k)<<" "<<fshname(k)<<endl;
    switch (fsh_sel_opt(k))
    {
      case 1 : // Selectivity coefficients 
      {
        if(phase_selcoff_fsh(k)<0)
        {
          write_input_log<<"Initial fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            log_selcoffs_fsh(k,jj)(1,nselages_in_fsh(k)) = log_selcoffs_fsh_in(k,jj)(1,nselages_in_fsh(k));
            write_input_log <<"Init coef:"<<endl<<exp(log_selcoffs_fsh(k,jj)(1,nselages_in_fsh(k))) <<endl;
          }
        }
      }
        break;
      case 2 : // Single logistic
      {
        if(phase_logist_fsh(k)<0)
        {
          logsel_slope_fsh(k,1) = logsel_slp_in_fsh(k,1)  ;
          write_input_log<<"Fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            logsel_slope_fsh(k,jj) = logsel_slp_in_fsh(k,jj)  ;
            sel50_fsh(k,jj)        =    sel_inf_in_fsh(k,jj)  ;
          }
        }
      }
      case 3 : // Double logistic 
      {
        if(phase_dlogist_fsh(k)<0)
        {
          write_input_log<<"Fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            logsel_slope_fsh(k,jj) = logsel_slp_in_fsh(k,jj)  ;
            sel50_fsh(k,jj)        =    sel_inf_in_fsh(k,jj)  ;
          }
        }
      }
      case 4 : // Selectivity spline initialize 
      /* {
        if(phase_sel_spl_fsh(k)<0)
        {
          write_input_log<<"Initial fishery spline to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            log_sel_spl_fsh(k,jj)(1,nnodes_tmp) = log_sel_spl_fsh_in(k,jj)(1,nnodes_tmp);
            // write_input_log <<"Init coef:"<<endl<<exp(log_sel_spl_fsh(k,jj)(1,nselages_in_fsh(k))) <<endl;
          }
          log_input(log_sel_spl_fsh);
        }
       }*/
     break;
    }
  }
  for (k=1;k<=nind;k++) 
  {
    write_input_log<<"Srvy sel phase: "<<phase_selcoff_ind(k)<<endl;
    if(phase_selcoff_ind(k)<0)
    {
      write_input_log<<"Fixing "<<indname(k)<<" indices sel to"<<endl<<n_sel_ch_ind(k)<<endl;
      for (int jj=1;jj<=n_sel_ch_ind(k);jj++) 
      {
        log_selcoffs_ind(k,jj)(1,nselages_in_ind(k)) = log_selcoffs_ind_in(k,jj)(1,nselages_in_ind(k));
        // write_input_log <<"Init coef:"<<endl<<exp(log_selcoffs_ind(k,jj)(1,nselages_in_ind(k))) <<endl;
      }
    }
    if(phase_logist_ind(k)<0)
    {
      write_input_log<<"Fixing index sel to"<<endl<<n_sel_ch_ind(k)<<endl;
      for (int jj=1;jj<=n_sel_ch_ind(k);jj++) 
      {
        logsel_slope_ind(k,jj) = logsel_slp_in_ind(k,jj)  ;
        // logsel_slope_ind(k,jj)    = 0.   ;
        sel50_ind(k,jj)           = sel_inf_in_ind(k,jj)  ;
      }
    }
  }
  log_input( logsel_slp_in_indv);
  write_input_log <<"Leaving parameter init secton"<<endl;
 END_CALCS

PRELIMINARY_CALCS_SECTION
  // Initialize age-specific changes in M if they are specified
  for (s=1;s<=nstk;s++)
  {
    dvar_vector yr_reg_st_tmp(1,nreg(s));
    yr_reg_st_tmp(1) = styr;
    for (r=2;r<=nreg(s);r++)
      yr_reg_st_tmp(r) = yy_shift_st(s,r);
    //Initialize matrix of M
    int r=1;
    for (i=styr;i<=endyr;i++)
    {
      if (i==yr_reg_st_tmp(r))
      {
        M(s,i) = Mest(mort_map(s,r));
        if (npars_Mage(mort_map(s,r))>0)
        {
          for (j=1;j<=npars_Mage(mort_map(s,r));j++)
            Mage_offset(mort_map(s,r),j) = Mage_offset_in(mort_map(s,r),j);
          int jj=1;
          for (j=1;j<=nages;j++)
          {
            if (j==ages_M_changes(mort_map(s,r),jj))
            {
              M(s,i,j) = natmortprior(mort_map(s,r))*mfexp(Mage_offset(mort_map(s,r),jj));
              jj++;
              if (npars_Mage(mort_map(s,r)) < jj) jj=npars_Mage(mort_map(s,r));
            }
            else
              if(j>1)
                M(s,i,j) = M(s,i,j-1);
          }
        }
        r++;
        if (nreg(s) < r) r=nreg(s);
      }
      else
        M(s,i) = M(s,i-1);
    }
  }
  log_input(M);
  Get_Age2length();
  if (use_vb_wt_age) // for now only uses 1st growth specification
  {
    for (k=1;k<=nfsh;k++)
      for ( i=styr; i<=endyr; i++ )
        wt_fsh(k,i) = value(wt_age_vb(1));
    for (k=1;k<=nind;k++)
      for ( i=styr; i<=endyr; i++ )
        wt_ind(k,i) = value(wt_age_vb(1));
    maturity(1) = value(maturity_vb(1));
  }

INITIALIZATION_SECTION
  Mest natmortprior; 
  steepness steepnessprior
  log_sigmar log_sigmarprior;

  log_Rzero    R_guess;
  mean_log_rec R_guess;
  
  log_Linf    log_Linfprior
  log_k       log_kprior
  log_Lo      log_Loprior
  log_sdage   log_sdageprior

  // log_avg_fmort -2.065
  log_q_ind log_qprior; 
  log_q_power_ind log_q_power_prior; 
  repl_F .1;

  sel50_fsh sel_inf_in_fshv 

  logsel_dslope_fsh logsel_dslp_in_fshv ;
  seld50_fsh sel_dinf_in_fshv 

  logsel_slope_ind logsel_slp_in_indv ;
  sel50_ind sel_inf_in_indv ;

  logsel_dslope_ind logsel_dslp_in_indv ;
  seld50_ind sel_dinf_in_indv ;

 //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==
PROCEDURE_SECTION
  fpen.initialize();
  for (k=1;k<=nind;k++) 
  {
    q_ind(k) = mfexp(log_q_ind(k) );
    q_power_ind(k) = mfexp(log_q_power_ind(k) );
  }

  // Main model calcs---------------------
  if(active(log_Linf(1))||active(log_k(1))||active(log_sdage(1))) //Ojo
    Get_Age2length();
  Get_Selectivity();
  Get_Mortality();
  Get_Bzero();
  Get_Numbers_at_Age();

  Get_Survey_Predictions();
  Get_Fishery_Predictions();
  // Objective function calcs------------
  evaluate_the_objective_function();
  if (last_phase())
    Get_Replacement_Yield();

  // Output calcs-------------------------
  if (sd_phase())
  {
    compute_spr_rates();
    Calc_Dependent_Vars();
    if (mcmcmode)
    {
      // Calc_Dependent_Vars();
      mcflag   = 0;
      mcmcmode = 0;
    }
    else
    {
      if (mcflag)
        Calc_Dependent_Vars();
    }
  }
  // Other calcs-------------------------
  if (mceval_phase())
  {
    if (oper_mod)
      Oper_Model();
    else
    {
      compute_spr_rates();
      write_mceval();
    }
  }
  if (do_fmort) Profile_F();
 //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==
FUNCTION write_mceval
  if (mcmcmode != 3)
    write_mceval_hdr();
  mcmcmode = 3;
  mceval<< model_name         << " "  ;
  mceval<< obj_fun            << " "  ;
  mceval<< obj_comps          << " "  ;
  // mceval<< rec_dev_future << " "  ;
  // mceval<<endl;
  get_msy();
  Future_projections();
  Calc_Dependent_Vars();
  mceval<<
  q_ind(1,1)  << " "<< 
  M(1,endyr,1)<< " "<< 
  steepness << " "<< 
  depletion << " "<< 
  MSY       << " "<< 
  MSYL      << " "<< 
  Fmsy      << " "<< 
  Fcur_Fmsy << " "<< 
  Bcur_Bmsy << " "<< 
  Bmsy      << " "<< 
  ABCBiom   << " "<< 
  F35       << " "<<
  F40       << " "<<
  F50       << " "<<
  SSB_fut(1,endyr_fut) << " "<< 
  SSB_fut(2,endyr_fut) << " "<< 
  SSB_fut(3,endyr_fut) << " "<< 
  SSB_fut(4,endyr_fut) << " "<< 
  SSB_fut(5,endyr_fut) << " "<< 
  catch_future(1,1,styr_fut)    << " "<<  
  catch_future(1,2,styr_fut)    << " "<<  
  catch_future(1,3,styr_fut)    << " "<<  
  catch_future(1,4,styr_fut)    << " "<<  endl;
//-----TRANSFORMATION FUNCION AGE->LENGTH--------------------------------------------------
FUNCTION Get_Age2length
 // This subroutine allows convert an age composition to length composition. For example: if there is a matrix C(1,nyears,1,nages), 
 // the vectorial operation: Cl=C*Prob_length,  returns a matrix Cl(1,nyears,1,nlength) whose sum over all the lengths is the same 
 // the sum over all age groups..
 // (Cristian Canales)
 // by default values
  // Linf=Linfprior;// Asymptotic length
  // k_coeff=kprior;
  // Lo=Loprior;// first length (corresponds to first age-group)
  // sdage=sdageprior;// coefficient of variation of length-at-age
 // if some of these are estimated.
  Linf    = mfexp(log_Linf);
  k_coeff = mfexp(log_k);
  Lo      = mfexp(log_Lo);
  sdage   = mfexp(log_sdage);
  for (r=1;r<=ngrowth;r++)
  {
    int i, j;
    mu_age(r,1)=Lo(r); // first length (modal)
    for (i=2;i<=nages;i++)
      mu_age(r,i) = Linf(r)*(1.-exp(-k_coeff(r)))+exp(-k_coeff(r))*mu_age(r,i-1); // the mean length by age group
    wt_age_vb(r) = lw_a * pow(mu_age(r), lw_b);
    maturity_vb(r) = 1/(1 + exp(32.93 - 1.45*mu_age(r)));
    sigma_age(r)=sdage(r)*mu_age(r); // standard deviation of length-at-age
    P_age2len(r) = ALK( mu_age(r), sigma_age(r), len_bins);
  }
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
  //RETURN_ARRAYS_INCREMENT();
  int i, j;
  dvariable z1;
  dvariable z2;
  int si,ni; si=mu.indexmin(); ni=mu.indexmax();
  int sj,nj; sj=x.indexmin(); nj=x.indexmax();
  dvar_matrix pdf(si,ni,sj,nj);
  double xs;
  pdf.initialize();
  for(i=si;i<=ni;i++) //loop over ages
  {
    for(j=sj;j<=nj;j++) //loop over length bins
    {
      if (j<nj)
        xs=0.5*(x[sj+1]-x[sj]);  // accounts for variable bin-widths...?
      z1=((x(j)-xs)-mu(i))/sig(i);
      z2=((x(j)+xs)-mu(i))/sig(i);
      pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
    }//end nbins
    pdf(i)/=sum(pdf(i));
  }//end nage
  //RETURN_ARRAYS_DECREMENT();
  return(pdf);
//---------------------------------------------------------------------------
FUNCTION Get_Replacement_Yield
  // compute next year's yield and SSB and add penalty to ensure F gives same SSB... 
  dvar_matrix ntmp(1,nstk,1,nages);
  for (s=1;s<=nstk;s++)
    ntmp(s) = natage(s,endyr+1);
  dvariable SSBnext;
  dvar_matrix Ftmp(1,nfsh,1,nages);
  dvar_matrix Ctmp(1,nstk,1,nages);
  dvar_matrix Ztmp(1,nstk,1,nages);
  dvar_matrix Stmp(1,nstk,1,nages);
  Ctmp.initialize();
  for (s=1;s<=nstk;s++)
    Ztmp(s)  = M(s,endyr);
  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
    sumF(sel_map(1,k)) += sum(F(k,endyr));
  for (k=1;k<=nfsh;k++)
  {
    int istk   = sel_map(1,k);
    Ftmp(k)    = repl_F(istk)*sum(F(k,endyr)) / sumF(istk);
    Ztmp(istk)+= Ftmp(k);
  }
  Stmp = mfexp(-Ztmp);
  for (k=1;k<=nfsh;k++)
  {
    int istk = sel_map(1,k);
    Ctmp(istk) += elem_prod(wt_fsh(k,endyr),elem_prod(elem_div(Ftmp(k),Ztmp(istk)),elem_prod(1.-Stmp(istk),ntmp(istk))) );
  }
  for (s=1;s<=nstk;s++)
  {
    repl_yld(s) = sum(Ctmp(s)) ;
    ntmp(s)(2,nages) = ++elem_prod(Stmp(s)(1,nages-1),ntmp(s)(1,nages-1));
    ntmp(s,nages)   += ntmp(s,nages)*Stmp(s,nages);
    ntmp(s,1)        = mean(mod_rec(s));
    repl_SSB(s)  = elem_prod(ntmp(s), pow(Stmp(s),spmo_frac)) * wt_mature(s); 
    obj_fun  += 200.*square(log(Sp_Biom(s,endyr))-log(repl_SSB(s)));
  }
  
FUNCTION Get_Selectivity
  // Calculate the logistic selectivity (Only if being used...)   
  for (k=1;k<=nfsh;k++)
  {
    switch (fsh_sel_opt(k))
    {
      case 1 : // Selectivity coefficients 
      //---Calculate the fishery selectivity from the sel_coffs (Only if being used...)   
      {
        int isel_ch_tmp = 1 ;
        dvar_vector sel_coffs_tmp(1,nselages_fsh(k));
        for (i=styr;i<=endyr;i++)
        {
          if (i==yrs_sel_ch_fsh(k,isel_ch_tmp)) 
          {
            sel_coffs_tmp.initialize();
            sel_coffs_tmp = log_selcoffs_fsh(k,isel_ch_tmp);
            avgsel_fsh(k,isel_ch_tmp)              = log(mean(mfexp(sel_coffs_tmp)));
            // Increment if there is still space to do so...
            if (isel_ch_tmp<n_sel_ch_fsh(k))
              isel_ch_tmp++;
          }
         // Need to flag for changing selectivity....XXX
          log_sel_fsh(k,i)(1,nselages_fsh(k))        = sel_coffs_tmp;
          log_sel_fsh(k,i)(nselages_fsh(k),nages)    = log_sel_fsh(k,i,nselages_fsh(k));
          log_sel_fsh(k,i)                                  -= log(mean(mfexp(log_sel_fsh(k,i) )));
        }
      }
      break;
      case 2 : // Single logistic
      {
        sel_slope_fsh(k) = mfexp(logsel_slope_fsh(k));
        int isel_ch_tmp = 1 ;
        dvariable sel_slope_tmp = sel_slope_fsh(k,isel_ch_tmp);
        dvariable sel50_tmp     = sel50_fsh(k,isel_ch_tmp);
        for (i=styr;i<=endyr;i++)
        {
          if (i==yrs_sel_ch_fsh(k,isel_ch_tmp)) 
          {
            sel_slope_tmp = sel_slope_fsh(k,isel_ch_tmp);
            sel50_tmp     =     sel50_fsh(k,isel_ch_tmp);
            if (isel_ch_tmp<n_sel_ch_fsh(k))
              isel_ch_tmp++;
          }
          log_sel_fsh(k,i)(1,nselages_fsh(k))     = -1.*log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                                                ( age_vector(1,nselages_fsh(k)) - sel50_tmp) ));
          log_sel_fsh(k,i)(nselages_fsh(k),nages) = log_sel_fsh(k,i,nselages_fsh(k));
        }
    }
    break;
    case 3 : // Double logistic
    {
      sel_slope_fsh(k)  = mfexp(logsel_slope_fsh(k));
      sel_dslope_fsh(k) = mfexp(logsel_dslope_fsh(k));
      int isel_ch_tmp = 1 ;
      dvariable sel_slope_tmp = sel_slope_fsh(k,isel_ch_tmp);
      dvariable sel50_tmp     = sel50_fsh(k,isel_ch_tmp);
      dvariable sel_dslope_tmp = sel_dslope_fsh(k,isel_ch_tmp);
      dvariable seld50_tmp     = seld50_fsh(k,isel_ch_tmp);
      for (i=styr;i<=endyr;i++)
      {
        if (i==yrs_sel_ch_fsh(k,isel_ch_tmp)) 
        {
          sel_slope_tmp  = sel_slope_fsh(k,isel_ch_tmp);
          sel50_tmp      =     sel50_fsh(k,isel_ch_tmp);
          sel_dslope_tmp = sel_dslope_fsh(k,isel_ch_tmp);
          seld50_tmp     =     seld50_fsh(k,isel_ch_tmp);
          if (isel_ch_tmp<n_sel_ch_fsh(k))
            isel_ch_tmp++;
        }
        log_sel_fsh(k,i)(1,nselages_fsh(k))     =
                     -log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                     ( age_vector(1,nselages_fsh(k)) - sel50_tmp) ))+
                     log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                     ( age_vector(1,nselages_fsh(k)) - seld50_tmp))) );

        log_sel_fsh(k,i)(nselages_fsh(k),nages) = 
                     log_sel_fsh(k,i,nselages_fsh(k));

        log_sel_fsh(k,i) -= max(log_sel_fsh(k,i));  
      }
    }
    break;
    //---Calculate the fishery selectivity from the sel_spl from nodes...
    case 4 : // Splines
     break;
    } // End of switch for fishery selectivity type
  } // End of fishery loop
  // Survey specific---
  for (k=1;k<=nind;k++)
  {
    switch (ind_sel_opt(k))
    {
      case 1 : // Selectivity coefficients
      //---Calculate the fishery selectivity from the sel_coffs (Only if being used...)   
      {
        int isel_ch_tmp = 1 ;
        dvar_vector sel_coffs_tmp(1,nselages_ind(k));
        for (i=styr;i<=endyr;i++)
        {
          if (i==yrs_sel_ch_ind(k,isel_ch_tmp)) 
          {
            sel_coffs_tmp.initialize();
            sel_coffs_tmp = log_selcoffs_ind(k,isel_ch_tmp);
            avgsel_ind(k,isel_ch_tmp)              = log(mean(mfexp(sel_coffs_tmp)));
            if (isel_ch_tmp<n_sel_ch_ind(k))
              isel_ch_tmp++;
          }
          log_sel_ind(k,i)(1,nselages_ind(k))        = sel_coffs_tmp;
          log_sel_ind(k,i)(nselages_ind(k),nages)    = log_sel_ind(k,i,nselages_ind(k));
          log_sel_ind(k,i)                                  -= log(mean(mfexp(log_sel_ind(k,i)(q_age_min(k),q_age_max(k))))); 
        }
      }
  
        break;
      case 2 : // Asymptotic logistic
        {
          sel_slope_ind(k) = mfexp(logsel_slope_ind(k));
          int isel_ch_tmp = 1 ;
          dvariable sel_slope_tmp = sel_slope_ind(k,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_ind(k,isel_ch_tmp);
          for (i=styr;i<=endyr;i++)
          {
            if (i==yrs_sel_ch_ind(k,isel_ch_tmp)) 
            {
              sel_slope_tmp = sel_slope_ind(k,isel_ch_tmp);
              sel50_tmp     =     sel50_ind(k,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_ind(k))
                isel_ch_tmp++;
            }
            log_sel_ind(k,i) = - log( 1.0 + mfexp(-sel_slope_tmp * ( age_vector - sel50_tmp) ));
            // log_sel_ind(k,i)                                  -= log(mean(mfexp(log_sel_ind(k,i)(q_age_min(k),q_age_max(k))))); 
          }
        }
        break;
      case 3 : // Double logistic
        {
          sel_slope_ind(k)  = mfexp(logsel_slope_ind(k));
          sel_dslope_ind(k) = mfexp(logsel_dslope_ind(k));
          int isel_ch_tmp = 1 ;
          dvariable sel_slope_tmp = sel_slope_ind(k,isel_ch_tmp);
          dvariable sel50_tmp     = sel50_ind(k,isel_ch_tmp);
          dvariable sel_dslope_tmp = sel_dslope_ind(k,isel_ch_tmp);
          dvariable seld50_tmp     = seld50_ind(k,isel_ch_tmp);
          for (i=styr;i<=endyr;i++)
          {
            if (i==yrs_sel_ch_ind(k,isel_ch_tmp)) 
            {
              sel_slope_tmp  = sel_slope_ind(k,isel_ch_tmp);
              sel50_tmp      =     sel50_ind(k,isel_ch_tmp);
              sel_dslope_tmp = sel_dslope_ind(k,isel_ch_tmp);
              seld50_tmp     =     seld50_ind(k,isel_ch_tmp);
              if (isel_ch_tmp<n_sel_ch_ind(k))
                isel_ch_tmp++;
            }
            log_sel_ind(k,i)(1,nselages_ind(k))     =
                         -log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                         ( age_vector(1,nselages_ind(k)) - sel50_tmp) ))+
                         log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                         ( age_vector(1,nselages_ind(k)) - seld50_tmp))) );

            log_sel_ind(k,i)(nselages_ind(k),nages) = 
                         log_sel_ind(k,i,nselages_ind(k));

            log_sel_ind(k,i) -= max(log_sel_ind(k,i));  
            log_sel_ind(k,i)                                  -= log(mean(mfexp(log_sel_ind(k,i)(q_age_min(k),q_age_max(k))))); 
          }
        }
      break;
    }// end of swtiches for indices selectivity
  } // End of indices loop

  // Map selectivities across fisheries and indices as needed.
  for (k=1;k<=nfsh;k++)
    if (sel_map(3,k)!=k)  // If 3rd row shows a different fishery then use that fishery
      log_sel_fsh(k) = log_sel_fsh(sel_map(3,k));

  for (k=1+nfsh;k<=nfsh_and_ind;k++)
    if (sel_map(2,k)!=2) 
      log_sel_ind(k-nfsh) = log_sel_fsh(sel_map(3,k));
    else if (sel_map(3,k)!=(k-nfsh)) 
      log_sel_ind(k-nfsh) = log_sel_ind(sel_map(3,k));

  sel_fsh = mfexp(log_sel_fsh);
  sel_ind = mfexp(log_sel_ind);

FUNCTION Get_NatMortality
  for (s=1;s<=nstk;s++)
  {
    dvar_vector yr_reg_st_tmp(1,nreg(s));
    yr_reg_st_tmp(1) = styr;
    for (r=2;r<=nreg(s);r++)
      yr_reg_st_tmp(r) = yy_shift_st(s,r);
    int r=1;
    int iyr_rw=1;
    for (i=styr;i<=endyr;i++)
    {
      if (i==yr_reg_st_tmp(r))
      {
        if (active(Mest(mort_map(s,r))))
        {
          natmort(s,i) = Mest(mort_map(s,r));
          M(s,i) = Mest(mort_map(s,r));
        }
        // Age varying part
        if (npars_Mage(mort_map(s,r))>0 && (active(Mest(mort_map(s,r))) || active(Mage_offset(mort_map(s,r)))))
        {
          int jj=1;
          for (j=1;j<=nages;j++)
          {
            if (j==ages_M_changes(mort_map(s,r),jj))
            {
              M(s,i,j) = natmortprior(mort_map(s,r))*mfexp(Mage_offset(mort_map(s,r),jj));
              jj++;
              if (npars_Mage(mort_map(s,r)) < jj) jj=npars_Mage(mort_map(s,r));
            }
            else
            {
              if(j>1)
                M(s,i,j) = M(s,i,j-1);
            }
          }
        }
        r++;
        if (nreg(s) < r) r=nreg(s);
      }
      else
      {
        natmort(s,i) = Mest(mort_map(s,yy_sr(s,i)));
        // Time varying part
        if (npars_rw_M(s)>0 && active(M_rw(s)))
        {
          if (i==yrs_rw_M(s,iyr_rw))
          {
            M(s,i) = M(s,i-1)*mfexp(M_rw(s,iyr_rw));
            iyr_rw++;
            if (npars_rw_M(s) < iyr_rw) 
										iyr_rw=npars_rw_M(s);
          }
          else
            M(s,i) = M(s,i-1);
        }
        else
          M(s,i) = M(s,i-1);
      }
    }
  }
	// natmort = M;

FUNCTION Get_Mortality2
  Get_NatMortality();
  Z       = M;
  for (k=1;k<=nfsh;k++)
  {
    F(k)   = elem_div(catage(k),natage(sel_map(1,k)));
    Z(sel_map(1,k))     += F(k);
  }
  S = mfexp(-1.*Z);

FUNCTION Get_Mortality
  Get_NatMortality();
  Z = M;
  if (!Popes)
  {
    Fmort.initialize();
    for (k=1;k<=nfsh;k++)
    {
      Fmort(sel_map(1,k)) +=  fmort(k);
      for (i=styr;i<=endyr;i++)
      {
        F(k,i)   =  fmort(k,i) * sel_fsh(k,i) ;
        Z(sel_map(1,k),i)    += F(k,i);
      }
    }
    S  = mfexp(-1.*Z);
  }
  

FUNCTION Get_Numbers_at_Age
  // natage(s,styr,1) = mfexp(mean_log_rec(cum_regs(s)+yy_sr(s,styr)) + rec_dev(s,styr)); 
  // Recruitment in subsequent years
  for (s=1;s<=nstk;s++)
  {
    for (i=styr+1;i<=endyr;i++)
      natage(s,i,1)=mfexp(mean_log_rec(cum_regs(s)+yy_sr(s,i))+rec_dev(s,i));
    
    mod_rec(s,styr)  = natage(s,styr,1);
    
    for (i=styr;i<=endyr;i++)
    {
      if (Popes)
      {
        dvariable  t1=mfexp(-natmort(s,i)*0.5);
        dvariable  t2=mfexp(-natmort(s,i));
        Catch_at_Age(s,i);
        // Pope's approximation //   Next year N     =   This year x NatSurvivl - catch
        natage(s,i+1)(2,nages) = ++(natage(s,i)(1,nages-1)*t2 - catage_tot(s,i)(1,nages-1)*t1);
        Ftot(s,i)(1,nages-1) = log(natage(s,i)(1,nages-1)) - --log(natage(s,i+1)(2,nages)) - natmort(s,i);
        natage(s,i+1,nages)   += natage(s,i,nages)*t2 - catage_tot(s,i,nages)*t1;
        // Approximation to "F" continuous form for computing within-year sp biomass
        Ftot(s,i,nages)      = log(natage(s,i,nages-1)+natage(s,i,nages)) -log(natage(s,i+1,nages)) -natmort(s,i);
        dvariable ctmp=sum(catage_tot(s,i));
        for (k=1;k<=nfsh;k++)
        {
          if (sel_map(1,k) == s)
            F(k,i)  = Ftot(s,i) * sum(catage(k,i))/ctmp;
        }
        Z(s,i)    = Ftot(s,i)+natmort(s,i);
        S(s,i)    = mfexp(-Z(s,i));
      }
      else // Baranov
      {
        // get_Fs( i ); //ojo, add switch here for different catch equation XX
        // if (i!=endyr)
        // {
          natage(s,i+1)(2,nages) = ++elem_prod(natage(s,i)(1,nages-1),S(s,i)(1,nages-1));
          natage(s,i+1,nages)   +=natage(s,i,nages)*S(s,i,nages);
        // }
      }
      Catch_at_Age(s,i);
      Sp_Biom(s,i)  = elem_prod(natage(s,i),pow(S(s,i),spmo_frac)) * wt_mature(s);
      if (i<endyr) mod_rec(s,i+1)  = natage(s,i+1,1);
    }
  }
  
  for (s=1;s<=nstk;s++)
  {
    for (r=2; r<=nreg(s); r++)
    {
      Bzero(cum_regs(s)+r) = Sp_Biom(s,reg_shift(s,r-1)-rec_age) ;  //reg_shift(s,i-1)-rec_age //(reg_shift(s,i-1)-nages)+1
    }
  }
  for (r=1; r<=nregs; r++)
  {
    phizero(r) = Bzero(r)/Rzero(r);

    int irec=rec_map(stk_reg_map(1,r),stk_reg_map(2,r));
    switch (SrType)
    {
      case 1:
        alpha(r) = log(-4.*steepness(irec)/(steepness(irec)-1.));
        break;
      case 2:
      {
        alpha(r)  =  Bzero(r) * (1. - (steepness(irec) - 0.2) / (0.8*steepness(irec)) ) / Rzero(r);
        beta(r)   = (5. * steepness(irec) - 1.) / (4. * steepness(irec) * Rzero(r));
      }
      break;
      case 4:
      {
        beta(r)  = log(5.*steepness(irec))/(0.8*Bzero(r)) ;
        alpha(r) = log(Rzero(r)/Bzero(r))+beta(r)*Bzero(r);
      }
        break;
    }
  }

FUNCTION Get_Survey_Predictions
  // Survey computations------------------
  dvariable sum_tmp;
  sum_tmp.initialize();
  int ii;
  int iyr;
  for (k=1;k<=nind;k++)
  {
    int istk = sel_map(1,k+nfsh);
    // Set rest of q's in time series equal to the random walk for current (avoids tricky tails...)
    for (i=2;i<=(1+npars_rw_q(k));i++)
    {
      // get index for the number of observations (can be different than number of q's)
      ii = yrs_rw_q(k,i-1) - yrs_ind(k,1) + 1;  
      q_ind(k,ii)  = q_ind(k,ii-1)*mfexp(log_rw_q_ind(k,i-1));
      for (iyr=ii+1;iyr<=nyrs_ind(k);iyr++)
        q_ind(k,iyr)  = q_ind(k,ii);
    }
    for (i=1;i<=nyrs_ind(k);i++)
    {        
      iyr=yrs_ind(k,i);
      pred_ind(k,i) = q_ind(k,i) * pow(elem_prod(natage(istk,iyr),pow(S(istk,iyr),ind_month_frac(k))) * 
                                     elem_prod(sel_ind(k,iyr) , wt_ind(k,iyr)),q_power_ind(k));
    }
    for (i=1;i<=nyrs_ind_age(k);i++)
    {        
      iyr = yrs_ind_age(k,i); 
      dvar_vector tmp_n   = elem_prod(pow(S(istk,iyr),ind_month_frac(k)),elem_prod(sel_ind(k,iyr),natage(istk,iyr)));  
      sum_tmp             = sum(tmp_n);
      if (use_age_err)
        eac_ind(k,i)      = age_err * tmp_n/sum_tmp;
      else
        eac_ind(k,i)      = tmp_n/sum_tmp;
    }
    dvar_vector tmp_n(1,nages);
    for (i=1;i<=nyrs_ind_length(k);i++)
    {        
      iyr          = yrs_ind_length(k,i); 
      tmp_n        = elem_prod(pow(S(istk,iyr),ind_month_frac(k)),elem_prod(sel_ind(k,iyr),natage(istk,iyr)));  
      sum_tmp      = sum(tmp_n);
      tmp_n       /= sum_tmp;
      int igrowth = growth_map(istk,yy_sr(istk,iyr));
      elc_ind(k,i) = tmp_n * P_age2len(igrowth) ;
    }
    iyr=yrs_ind(k,nyrs_ind(k));
    dvar_vector natagetmp = elem_prod(S(istk,endyr),natage(istk,endyr));
    natagetmp(2,nages) = ++natagetmp(1,nages-1)*1.;
    natagetmp(1)       = SRecruit(Sp_Biom(istk,endyr+1-rec_age),cum_regs(istk)+yy_sr(istk,endyr+1));
    natagetmp(nages)  += natage(istk,endyr,nages)*S(istk,endyr,nages);
    // Assume same survival in 1st part of next year as same as first part of current
    pred_ind_nextyr(k) = q_ind(k,nyrs_ind(k)) * pow(elem_prod(natagetmp,pow(S(istk,endyr),ind_month_frac(k))) * 
                                     elem_prod(sel_ind(k,endyr) , wt_ind(k,endyr)),q_power_ind(k));
  }

FUNCTION Get_Fishery_Predictions
  for (k=1;k<=nfsh;k++)
  {
    for (i=1; i<=nyrs_fsh_age(k); i++)
    {
      if (use_age_err)
        eac_fsh(k,i) = age_err * catage(k,yrs_fsh_age(k,i))/sum(catage(k,yrs_fsh_age(k,i)));
      else
        eac_fsh(k,i) = catage(k,yrs_fsh_age(k,i))/sum(catage(k,yrs_fsh_age(k,i)));
      eac_fsh(k,i) /= sum(eac_fsh(k,i));
    }

 // predicted length compositions !!
    for (i=1; i<=nyrs_fsh_length(k); i++)
    {
      int istk = sel_map(1,k);
      int iyr = yrs_fsh_length(k,i);
      int igrowth = growth_map(istk,yy_sr(istk,iyr));
      elc_fsh(k,i) = catage(k,yrs_fsh_length(k,i))*P_age2len(igrowth);
      elc_fsh(k,i) /= sum(elc_fsh(k,i));
    }
  }

FUNCTION Calc_Dependent_Vars
  get_msy();

  if (phase_proj>0) Future_projections();
  N_NoFsh.initialize();
  dvar_matrix Nnext(1,nstk,1,nages);
  Nnext.initialize();
  for (s=1;s<=nstk;s++)
  {
    N_NoFsh(s,styr) = natage(s,styr);
    for (i=styr_sp;i<=styr;i++)
      Sp_Biom_NoFish(s,i) = Sp_Biom(s,i);
    for (i=styr;i<=endyr;i++)
    {                 
      recruits(s,i)  = natage(s,i,1);
      if (i>styr)
      {
        N_NoFsh(s,i,1)        = recruits(s,i);
        N_NoFsh(s,i,1)       *= SRecruit(Sp_Biom_NoFish(s,i-rec_age),cum_regs(s)+yy_sr(s,i)) / SRecruit(Sp_Biom(s,i-rec_age),cum_regs(s)+yy_sr(s,i));
        N_NoFsh(s,i)(2,nages) = ++elem_prod(N_NoFsh(s,i-1)(1,nages-1),exp(-M(s,i-1)(1,nages-1)));
        N_NoFsh(s,i,nages)   += N_NoFsh(s,i-1,nages)*exp(-M(s,i-1,nages));
      }
      totbiom_NoFish(s,i) = N_NoFsh(s,i)*wt_pop(s);
      totbiom(s,i)        = natage(s,i)*wt_pop(s);
      Sp_Biom_NoFish(s,i) = N_NoFsh(s,i)*elem_prod(pow(exp(-M(s,i)),spmo_frac) , wt_mature(s)); 
      Sp_Biom_NoFishRatio(s,i) = Sp_Biom(s,i) / Sp_Biom_NoFish(s,i) ;
      depletion(s)         = totbiom(s,endyr)/totbiom(s,styr);
      depletion_dyn(s)     = totbiom(s,endyr)/totbiom_NoFish(s,endyr);
    }
    B100(s) = phizero(cum_regs(s)+yy_sr(s,styr)) * mean(recruits(s)(styr_rec_est(s,1),endyr_rec_est(s,nreg(s)))); //Ojo
    //dvar_vector Nnext(1,nages);
    Nnext(s)(2,nages) = ++elem_prod(natage(s,endyr)(1,nages-1),S(s,endyr)(1,nages-1));
    Nnext(s,nages)   += natage(s,endyr,nages)*S(s,endyr,nages);
    // Compute SSB in next year using mean recruits for age 1 and same survival as in endyr
    Nnext(s,1)        = mfexp(mean_log_rec(cum_regs(s)+yy_sr(s,endyr+1))+rec_dev_future(s,endyr+1));
    Sp_Biom(s,endyr+1)= elem_prod(Nnext(s),pow(S(s,endyr),spmo_frac)) * wt_mature(s); 
    // Nnext(1)       = SRecruit(Sp_Biom(endyr+1-rec_age));
    ABCBiom(s)       = Nnext(s)*wt_pop(s);
    recruits(s,endyr+1) = Nnext(s,1);
    totbiom(s,endyr+1)  = ABCBiom(s);
  }
  // Now do OFL for next year...
  dvar_matrix seltmp(1,nfsh,1,nages);
  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_matrix Ztmp(1,nstk,1,nages);
  seltmp.initialize();
  Fatmp.initialize();
  Ztmp.initialize();
  for (k=1;k<=nfsh;k++)
    seltmp(k) = (sel_fsh(k,endyr));
  for (s=1;s<=nstk;s++)
    Ztmp(s) = (M(s,styr));
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k)            = (Fratio(k) * Fmsy(sel_map(1,k)) * seltmp(k));
    Ztmp(sel_map(1,k)) += Fatmp(k);
  }
  dvar_matrix survmsy(1,nstk,1,nages);
  survmsy = exp(-Ztmp);
  dvar_vector ctmp(1,nages);
  ctmp.initialize();
  OFL.initialize();
  for (k=1;k<=nfsh;k++)
  {
    for ( j=1; j<=nages; j++ )
      ctmp(j) = Nnext(sel_map(1,k),j) * Fatmp(k,j) * (1. - survmsy(sel_map(1,k),j)) / Ztmp(sel_map(1,k),j);
    OFL(sel_map(1,k)) += wt_fsh(k,endyr) * ctmp;
  }
  
FUNCTION void Catch_at_Age(const int& s, const int& i)
  dvariable vbio=0.;
  dvariable pentmp;
  dvar_vector Nmid(1,nages);
  dvar_vector Ctmp(1,nages);
  catage_tot(s,i).initialize();
  if (Popes)
  {
    Nmid = elem_prod(natage(s,i),mfexp(-M(s,i)/2) ); 
  }
  for (k=1;k<=nfsh;k++)
  {
    if (sel_map(1,k) == s)
    {
      if (Popes)
      {
        pentmp=0.;
        Ctmp = elem_prod(Nmid,sel_fsh(k,i));
        vbio = Ctmp*wt_fsh(k,i);
        //Kludge to go here...
        // dvariable SK = posfun( (.98*vbio - catch_bio(k,i))/vbio , 0.1 , pentmp );
        dvariable SK = posfun( (vbio - catch_bio(k,i))/vbio , 0.1 , pentmp );
        catch_tmp    = vbio - SK*vbio; 
        hrate        = catch_tmp / vbio;
        fpen(4) += pentmp;
        Ctmp *= hrate;                          
        if (hrate>1) {cout << catch_tmp<<" "<<vbio<<endl;exit(1);}
        catage_tot(s,i) += Ctmp;                      
        catage(k,i)    = Ctmp;                      
        if (last_phase())
          pred_catch(k,i) = Ctmp*wt_fsh(k,i);
      }
      else
      {
        catage(k,i) = elem_prod(elem_div(F(k,i),Z(s,i)),elem_prod(1.-S(s,i),natage(s,i)));
        pred_catch(k,i) = catage(k,i)*wt_fsh(k,i);
      }
    }
  }
  //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==
FUNCTION evaluate_the_objective_function
  // if (active(fmort_dev))   
  if (active(fmort))   
  {
    Cat_Like();
    Fmort_Pen();
  }
  Rec_Like();
  if (active(rec_dev))
    Age_Like();
  Srv_Like();
  Sel_Like();
  Compute_priors();
  for (r=1;r<=nregs;r++)
  {
    if (active(log_Rzero(r)))
      obj_fun += .5 * square(log_Rzero(r)-mean_log_rec(r));
  }
  //obj_fun += sum(.5 * square(log_Rzero(1)-mean_log_rec(1))); // A slight penalty to keep Rzero in reality...
  
  obj_comps.initialize();
  obj_comps(1)  = sum(catch_like);
  obj_comps(2)  = sum(age_like_fsh);
//------------------------------------------NEW-------------
  obj_comps(3)  = sum(length_like_fsh);
//-----------------------------------------------------------
  obj_comps(4)  = sum(sel_like_fsh);
  obj_comps(5)  = sum(ind_like);
  obj_comps(6)  = sum(age_like_ind);
  obj_comps(7)  = sum(length_like_ind);
  obj_comps(8)  = sum(sel_like_ind);
  obj_comps(9)  = sum(rec_like);
  obj_comps(10) = sum(fpen);
  obj_comps(11) = sum(post_priors_indq);
  obj_comps(12) = sum(post_priors);
  obj_fun     += sum(obj_comps);
FUNCTION Cat_Like
  // Eases into the catch-biomass likelihoods.  If too far off to start, full constraint to fit can be too aggressive
  catch_like.initialize();
  dvariable catch_pen;
  switch (current_phase())
  {
    case 1:
      catch_pen = .1;
      break;
    case 2:
      catch_pen = .5;
      break;
    case 3:
      catch_pen = .8;
      break;
    case 4:
      catch_pen = 1.0;
      break;
    case 5:
      catch_pen = 1;
      break;
    default:
      catch_pen = 1;
      break;
  }
  if (current_phase()>3)
  {
    for (k=1;k<=nfsh;k++)
      for (i=styr;i<=endyr;i++)
         catch_like(k) += .5*square(log(catch_bio(k,i)+.0001) - log(pred_catch(k,i)+.0001) )/catch_bio_lva(k,i);
  }
  else
  {
    for (k=1;k<=nfsh;k++)
      catch_like(k) += catchbiomass_pen * norm2(log(catch_bio(k)   
                      +.000001) - log(pred_catch(k) +.000001));
  }

  catch_like *= catch_pen;
FUNCTION Rec_Like
  rec_like.initialize();
  pred_rec.initialize();
  if (active(rec_dev))
  {
    sigmar = mfexp(log_sigmar);
    for (r=1;r<=nregs;r++)
    {
      int istk    = stk_reg_map(1,r);
      int ireg    = stk_reg_map(2,r);
      sigmarsq(r) = square(sigmar(rec_map(istk,ireg)));
    }
    if (current_phase()>2)
    {
      for (s=1;s<=nstk;s++)
      {
        for (r=1;r<=nreg(s);r++)
        {
          if (last_phase())
            pred_rec(s)(yy_shift_st(s,r),yy_shift_end(s,r)) = SRecruit(Sp_Biom(s)(yy_shift_st(s,r)-rec_age,yy_shift_end(s,r)-rec_age).shift(yy_shift_st(s,r))(yy_shift_st(s,r),yy_shift_end(s,r)),cum_regs(s)+r)(yy_shift_st(s,r),yy_shift_end(s,r));
          else 
            pred_rec(s)(yy_shift_st(s,r),yy_shift_end(s,r)) = .1+SRecruit(Sp_Biom(s)(yy_shift_st(s,r)-rec_age,yy_shift_end(s,r)-rec_age).shift(yy_shift_st(s,r))(yy_shift_st(s,r),yy_shift_end(s,r)),cum_regs(s)+r)(yy_shift_st(s,r),yy_shift_end(s,r));
        }
      }

      dvar_vector SSQRec(1,nregs);
      SSQRec.initialize();
      dvar_matrix chi(1,nstk,styr,endyr);
      chi.initialize();
      for (r=1;r<=nregs;r++)
      {
        int istk = stk_reg_map(1,r);
        int ireg = stk_reg_map(2,r);
        for (j=1;j<=nrecs_est_shift(r);j++)
        {
          chi(istk,yr_rec_est(r,j)) = log(mod_rec(istk,yr_rec_est(r,j))) - log(pred_rec(istk,yr_rec_est(r,j)));
        }
        SSQRec(r)   = norm2( chi(istk)(styr_rec_est(istk,ireg),endyr_rec_est(istk,ireg)) ) ;
        m_sigmarsq(r) =  SSQRec(r)/nrecs_est_shift(r);
        m_sigmar(r)   =  sqrt(m_sigmarsq(r));
      }

      if (current_phase()>4||last_phase())
        for (r=1;r<=nregs;r++)
        {
          int istk = stk_reg_map(1,r);
          int ireg = stk_reg_map(2,r);
          rec_like(istk,1) += (SSQRec(r)+ m_sigmarsq(r)/2.)/(2*sigmarsq(r)) + nrecs_est_shift(r)*log_sigmar(rec_map(istk,ireg)); 
        }
      else
        for (r=1;r<=nregs;r++)
        {
          int istk = stk_reg_map(1,r);
          int ireg = stk_reg_map(2,r);
          rec_like(istk,1) += .1*((SSQRec(r)+ m_sigmarsq(r)/2.)/(2*sigmarsq(r)) + nrecs_est_shift(r)*log_sigmar(rec_map(istk,ireg))); 
        }
    }

    if (last_phase())
    {
      for (s=1;s<=nstk;s++)
      {
        // Variance term for the parts not estimated by sr curve
        if ( styr_rec_est(s,1) > styr_rec )
          rec_like(s,4) += .5*norm2( rec_dev(s)(styr_rec,styr_rec_est(s,1)-1) )/sigmarsq(cum_regs(s)+1) + ((styr_rec_est(s,1)-1)-styr_rec)*log(sigmar(rec_map(s,1))) ;
        
        if ( endyr > endyr_rec_est(s,nreg(s)) )
          rec_like(s,4) += .5*norm2( rec_dev(s)(endyr_rec_est(s,nreg(s))+1,endyr) )/sigmarsq(cum_regs(s)+nreg(s)) + (endyr-(endyr_rec_est(s,nreg(s))+1))*log(sigmar(rec_map(s,nreg(s)))) ;
        
        for (r=2;r<=nreg(s);r++)
        {
          if ( (styr_rec_est(s,r)-1) > endyr_rec_est(s,r-1) )
            rec_like(s,4) += .5*norm2( rec_dev(s)(endyr_rec_est(s,r-1)+1,styr_rec_est(s,r)-1) )/sigmarsq(cum_regs(s)+(r-1)) + ((styr_rec_est(s,r)-1)-(endyr_rec_est(s,r-1)+1))*log(sigmar(rec_map(s,r-1))) ; // Ojo
        }
        
        for (r=1;r<=nreg(s);r++)
        {
          int iregs = cum_regs(s)+r;
          for (j=1;j<=(nrecs_est_shift(iregs)-1);j++)
          {
            if ((yr_rec_est(iregs,j+1)-yr_rec_est(iregs,j)) > 1)
              rec_like(s,4) += .5*norm2( rec_dev(s)(yr_rec_est(iregs,j)+1,yr_rec_est(iregs,j+1)-1) )/sigmarsq(iregs) + ((yr_rec_est(iregs,j+1)-1)-(yr_rec_est(iregs,j)+1))*log(sigmar(rec_map(s,r))) ;
          }
        }
      }
    }
    else // JNI comment next line
    {
      for (r=1;r<=nregs;r++)
      {
        int istk = stk_reg_map(1,r);
        for (j=1;j<=nrecs_est_shift(r);j++)
        {
          rec_like(istk,2) += square(rec_dev(istk,yr_rec_est(r,j)));
        }
      }
    }

    for (s=1;s<=nstk;s++)
    {
      rec_like(s,2) += norm2( rec_dev(s)(styr_rec_est(s,1),endyr) ) ;
    }

    if (active(rec_dev_future))
    {
      // Future recruitment variability (based on past)
      for (s=1;s<=nstk;s++)
      {
        sigmar_fut(s)   = sigmar(rec_map(s,nreg(s))) ;
        rec_like(s,3) += norm2(rec_dev_future(s))/(2*square(sigmar_fut(s)))+ size_count(rec_dev_future(s))*log(sigmar_fut(s));
      }
    }
  }
FUNCTION Compute_priors
  post_priors.initialize();
  post_priors_indq.initialize();
  for (k=1;k<=nind;k++)
  {
    if (active(log_q_ind(k)))
      post_priors_indq(k) += square(log(q_ind(k,1)/qprior(k)))/(2.*cvqprior(k)*cvqprior(k)); 
    if (active(log_q_power_ind(k)))
      post_priors_indq(k) += square(log(q_power_ind(k)/q_power_prior(k)))/(2.*cvq_power_prior(k)*cvq_power_prior(k)); 
    if (active(log_rw_q_ind(k)))
      for (int i=1;i<=npars_rw_q(k);i++)
      {
        post_priors_indq(k) += square(log_rw_q_ind(k,i))/ (2.*sigma_rw_q(k,i)*sigma_rw_q(k,i)) ;
      }
     //  -q_power_prior(k))/(2*cvq_power_prior(k)*cvq_power_prior(k)); 
  }

  for (r=1;r<=nmort;r++)
    if (active(Mest(r)))
      post_priors(r,1) += square(log(Mest(r)/natmortprior(r)))/(2.*cvnatmortprior(r)*cvnatmortprior(r)); 

  for (r=1;r<=nmort;r++)
    if (active(Mage_offset(r)))
      post_priors(r,1) += norm2(Mage_offset(r))/(2.*cvnatmortprior(r)*cvnatmortprior(r)); 

  for (s=1;s<=nstk;s++)
    if (active(M_rw(s)))
      for (int i=1;i<=npars_rw_M(s);i++)
        post_priors(s,1) +=  square(M_rw(s,i))/ (2.*sigma_rw_M(s,i)*sigma_rw_M(s,i)) ;

  for (int r=1;r<=nrec;r++)
  {
    if (active(steepness(r)))
      post_priors(r,2) += square(log(steepness(r)/steepnessprior(r)))/(2*cvsteepnessprior(r)*cvsteepnessprior(r)); 
  }

  for (int r=1;r<=nrec;r++)
  {
    if (active(log_sigmar(r)))
      post_priors(r,3) += square(log(sigmar(r)/sigmarprior(r)))/(2*cvsigmarprior(r)*cvsigmarprior(r)); 
  }
//--------------------------NEW------------------------------------
  for (int r=1;r<=ngrowth;r++)
    if (active(log_Linf(r)))  
      post_priors(r,4) += square(log_Linf(r)-log_Linfprior(r))/(2*cvLinfprior(r)*cvLinfprior(r)); 

  for (int r=1;r<=ngrowth;r++)
    if (active(log_k(r)))
      post_priors(r,5) += square(log_k(r)-log_kprior(r))/(2*cvkprior(r)*cvkprior(r)); 

  for (int r=1;r<=ngrowth;r++)
    if (active(log_Lo(r)))
      post_priors(r,6) += square(log_Lo(r)-log_Loprior(r))/(2*cvLoprior(r)*cvLoprior(r)); 

  for (int r=1;r<=ngrowth;r++)
    if (active(log_sdage(r)))
      post_priors(r,7) += square(log_sdage(r)-log_sdageprior(r))/(2*cvsdageprior(r)*cvsdageprior(r)); 
FUNCTION Fmort_Pen
  // Phases less than 3, penalize High F's---------------------------------
  if (current_phase()<3)
    fpen(1) += 1.* norm2(F - .2);
  else 
    fpen(1) += 0.0001*norm2(F - .2); 

  // for (k=1;k<=nfsh;k++)  fpen(2) += 20.*square(mean(fmort_dev(k)) ); // this is just a normalizing constraint (fmort_devs sum to zero) }
FUNCTION Sel_Like 
  sel_like_fsh.initialize();
  sel_like_ind.initialize();
  for (k=1;k<=nfsh;k++)
  {
    if (active(logsel_slope_fsh(k)))
    {
      for (i=2;i<=n_sel_ch_fsh(k);i++)
      {
          int iyr = yrs_sel_ch_fsh(k,i) ;
          dvariable var_tmp = square(sel_sigma_fsh(k,i));
          sel_like_fsh(k,2)    += .5*norm2( log_sel_fsh(k,iyr-1) - log_sel_fsh(k,iyr) ) / var_tmp ;
      }
    }

    if (active(log_selcoffs_fsh(k)))
    {
      for (i=1;i<=n_sel_ch_fsh(k);i++)
      {
        int iyr = yrs_sel_ch_fsh(k,i) ;
        // If curvature penalty is assumed....
        sel_like_fsh(k,1) += curv_pen_fsh(k)*norm2(first_difference( first_difference(log_sel_fsh(k,iyr))));
        // If curvature penalty (sigma) is estimated....
        // dvariable var=mfexp(2.0*logSdsmu_fsh(k));
        // sel_like_fsh(k,1) += 0.5*(size.count(log_sel_fsh(k,iyr))*log(var) +  norm2(first_difference( first_difference(log_sel_fsh(k,iyr)))) /var);
        if (i>1)
        {
          // This part is the penalty on the change itself--------------
          dvariable var_tmp = square(sel_sigma_fsh(k,i));
          sel_like_fsh(k,2)    += .5*norm2( log_sel_fsh(k,iyr-1) - log_sel_fsh(k,iyr) ) / var_tmp ;
        }
        int nagestmp = nselages_fsh(k);
        for (j=seldecage;j<=nagestmp;j++)
        {
          dvariable difftmp = log_sel_fsh(k,iyr,j-1)-log_sel_fsh(k,iyr,j) ;
          if (difftmp > 0.)
            sel_like_fsh(k,3)    += .5*square( difftmp ) / seldec_pen_fsh(k);
        }
        obj_fun            += 20 * square(avgsel_fsh(k,i)); // To normalize selectivities
      }
    }
  }
  for (k=1;k<=nind;k++)
  {
    if (active(logsel_slope_ind(k)))
    {
      for (i=2;i<=n_sel_ch_ind(k);i++)
      {
          int iyr = yrs_sel_ch_ind(k,i) ;
          dvariable var_tmp = square(sel_sigma_ind(k,i));
          sel_like_ind(k,2)    += .5*norm2( log_sel_ind(k,iyr-1) - log_sel_ind(k,iyr) ) / var_tmp ;
      }
    }
    if (active(log_selcoffs_ind(k)))
    {
      int nagestmp = nselages_ind(k);
      for (i=1;i<=n_sel_ch_ind(k);i++)
      {
        int iyr = yrs_sel_ch_ind(k,i) ;
        sel_like_ind(k,1) += curv_pen_ind(k)*norm2(first_difference( first_difference(log_sel_ind(k,iyr))));
        // This part is the penalty on the change itself--------------
        if (i>1)
        {
          dvariable var_tmp = square(sel_sigma_ind(k,i));
          sel_like_ind(k,2)    += .5*norm2( log_sel_ind(k,iyr-1) - log_sel_ind(k,iyr) ) / var_tmp ;
        }
        for (j=seldecage;j<=nagestmp;j++)
        {
          dvariable difftmp = log_sel_ind(k,iyr,j-1)-log_sel_ind(k,iyr,j) ;
          if (difftmp > 0.)
            sel_like_ind(k,3)    += .5*square( difftmp ) / seldec_pen_ind(k);
        }
        obj_fun            += 20. * square(avgsel_ind(k,i));  // To normalize selectivities
      }
    }
  }
FUNCTION Srv_Like
  // Fit to indices (log-Normal) -------------------------------------------
  ind_like.initialize();
  int iyr;
  for (k=1;k<=nind;k++)
    for (i=1;i<=nyrs_ind(k);i++)
    {
      // iyr = int(yrs_ind(k,i));
      ind_like(k) += square(log(obs_ind(k,i)) - log(pred_ind(k,i)) ) / 
                                   (2.*obs_lse_ind(k,i)*obs_lse_ind(k,i));
    }
  /* normal distribution option to add someday...
    for (i=1;i<=nyrs_ind(k);i++)
      ind_like(k) += square(obs_ind(k,i) - pred_ind(k,yrs_ind(k,i)) ) / 
                                   (2.*obs_se_ind(k,i)*obs_se_ind(k,i));
  */
FUNCTION Age_Like
  age_like_fsh.initialize();
  for (k=1;k<=nfsh;k++)
    for (int i=1;i<=nyrs_fsh_age(k);i++)
      age_like_fsh(k) -= n_sample_fsh_age(k,i)*(oac_fsh(k,i) + 0.001) * log(eac_fsh(k,i) + 0.001 ) ;
  age_like_fsh -= offset_fsh;
//-----------------------------------NEW-----------------------
  length_like_fsh.initialize();
  for (k=1;k<=nfsh;k++)
    for (int i=1;i<=nyrs_fsh_length(k);i++)
      length_like_fsh(k) -= n_sample_fsh_length(k,i)*(olc_fsh(k,i) + 0.001) * log(elc_fsh(k,i) + 0.001 ) ;
  length_like_fsh -= offset_lfsh;
//----------------------------------------------------------
  length_like_ind.initialize();
  for (k=1;k<=nind;k++)
    for (int i=1;i<=nyrs_ind_length(k);i++)
      length_like_ind(k) -= n_sample_ind_length(k,i)*(olc_ind(k,i) + 0.001) * log(elc_ind(k,i) + 0.001 ) ;
  length_like_ind -= offset_lind;
//----------------------------------------------------------
  age_like_ind.initialize();
  for (k=1;k<=nind;k++)
    for (int i=1;i<=nyrs_ind_age(k);i++)
      age_like_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(eac_ind(k,i) + 0.001 ) ;
  age_like_ind -= offset_ind;

FUNCTION Oper_Model
 // Initialize things used here only
  mc_count++;
  get_msy();
  Write_SimDatafile();
  Write_Datafile();
  dmatrix new_ind(1,nind,1,nyrs_ind);
  new_ind.initialize();

  int nsims;
  ifstream sim_in("nsims.dat");
  sim_in >> nsims; sim_in.close();

  dvector ran_ind_vect(1,nind);
  ofstream SaveOM("Om_Out.dat",ios::app);
  dvector C_tmp(1,nstk);
  dvector Fnow(1,nstk);
  for (s=1;s<=nstk;s++)
  {
    // Initialize recruitment in first year
    for (i=styr_fut-rec_age;i<styr_fut;i++)
      Sp_Biom_future(s,i) = Sp_Biom(s,i);
    nage_future(s,styr_fut)(2,nages)              = ++elem_prod(natage(s,endyr)(1,nages-1),S(s,endyr)(1,nages-1));
    nage_future(s,styr_fut,nages)                += natage(s,endyr,nages)*S(s,endyr,nages);

    // assume survival same as in last year...
    Sp_Biom_future(s,styr_fut) = elem_prod(nage_future(s,styr_fut),pow(S(s,endyr),spmo_frac)) * wt_mature(s); 
  }
  for (int isim=1;isim<=nsims;isim++)
  {
    cout<<isim<<" "<<cmp_no<<" "<<mc_count<<" "<<endl;
    // Copy file to get mean for Mgt Strategies
    system("init_stuff.bat");
    for (i=styr_fut;i<=endyr_fut;i++)
    {
      // Some unit normals...for generating data
      ran_ind_vect.fill_randn(rng);
      cout<<ran_ind_vect<<endl;
      // Create new indices observations
      // for (k = 1 ; k<= nind ; k++) new_ind(k) = mfexp(ran_ind_vect(k)*.2)*value(nage_future(i)*q_ind(k,nyrs_ind(k))*sel_ind(k,endyr)); // use value function since converts to a double
      // new_ind(1) = mfexp(ran_ind_vect(1)*0.2)*value(sum(nage_future(i)*q_ind(1,nyrs_ind(1))));
      if(styr_fut==i)
        new_ind(1) = mfexp(ran_ind_vect(1)*0.2)*value(wt_ind(1,endyr)*(natage(1,i-1))); //Ojo
      else
        new_ind(1) = mfexp(ran_ind_vect(1)*0.2)*value(wt_ind(1,endyr)*(nage_future(1,i-1))); //Ojo
      // now for Selecting which MP to use
      // Append new indices observation to datafile
      ifstream tacin("ctac.dat");
      int nobstmp;
      tacin >> nobstmp ;
      dvector t_tmp(1,nobstmp);
      tacin >> t_tmp;
      tacin.close();
      ofstream octac("ctac.dat");
      octac<<nobstmp+1<<endl;
      octac<<t_tmp<<endl;
      octac<<new_ind(1)<<endl;
      octac.close();
      system("ComputeTAC.bat " + (itoa(cmp_no,10))); // commandline function to get TAC (catchnext.dat)
     // Now read in TAC (actual catch)
     ifstream CatchNext("CatchNext.dat");
     CatchNext >> C_tmp; 
     CatchNext.close();
     //if (cmp_no==5) C_tmp=value((natmort(styr))*mean(t_tmp(nobstmp-2,nobstmp)));
     //if (cmp_no==6) C_tmp=value((natmort(styr))*.75*mean(t_tmp(nobstmp-2,nobstmp)));
     if (cmp_no==5) 
     {
       for (s=1;s<=nstk;s++)
         C_tmp(s) = min(C_tmp(s)*1.1,value((natmort(s,styr)*t_tmp(nobstmp))));
       ofstream cnext("CatchNext.dat");
       cnext <<C_tmp<<endl;
       cnext.close();
     }
     if (cmp_no==6) 
     {
       for (s=1;s<=nstk;s++)
         C_tmp(s) = min(C_tmp(s)*1.1,value(natmort(s,styr)*.75*t_tmp(nobstmp)));
       ofstream cnext("CatchNext.dat");
       cnext <<C_tmp<<endl;
       cnext.close();
     }

     for (s=1;s<=nstk;s++)
       Fnow(s) = value(SolveF2(endyr,nage_future(s,i), C_tmp(s),s));

      F_future(1,i) = sel_fsh(1,endyr) * Fnow; //Ojo
      for (s=1;s<=nstk;s++)
      {
        //Z_future(i)   = F_future(1,i) + max(natmort);
        Z_future(s,i)   = F_future(1,i) + mean(M(s));
        S_future(s,i)   = mfexp(-Z_future(s,i));
        nage_future(s,i,1)  = SRecruit( Sp_Biom_future(s,i-rec_age),cum_regs(s)+yy_sr(s,i) ) * mfexp(rec_dev_future(s,i)) ;     
        Sp_Biom_future(s,i) = wt_mature(s) * elem_prod(nage_future(s,i),pow(S_future(s,i),spmo_frac)) ;
        // Now graduate for the next year....
        if (i<endyr_fut)
        {
          nage_future(s,i+1)(2,nages) = ++elem_prod(nage_future(s,i)(1,nages-1),S_future(s,i)(1,nages-1));
          nage_future(s,i+1,nages)   += nage_future(s,i,nages)*S_future(s,i,nages);
        }
        catage_future(s,i) = 0.; 
        for (k = 1 ; k<= nfsh ; k++)
          if (sel_map(1,k) == s)
            catage_future(s,i) += elem_prod(nage_future(s,i) , elem_prod(F_future(k,i) , elem_div( ( 1.- S_future(s,i) ) , Z_future(s,i))));
      }
  
      for (s=1;s<=nstk;s++)
        SaveOM << model_name       <<
          " "  << cmp_no           <<
          " "  << mc_count         <<
          " "  << isim             <<
          " "  << i                <<
          " "  << s                <<
          " "  << Fnow(s)          <<
          " "  << Fnow(s)/Fmsy(s)  <<
          " "  << Sp_Biom_future(s,i-rec_age)         <<
          " "  << nage_future(s,i)                    <<
          " "  << catage_future(s,i)*wt_fsh(1,endyr)  <<
          " "  << mean(M(s))                          <<
          " "  << t_tmp(nobstmp)                      <<
        endl;
    }
  }
  // if (mc_count>5) exit(1);
  SaveOM.close();
  if (!mceval_phase())
    exit(1);

FUNCTION void get_future_Fs(const int& s,const int& i,const int& iscenario)
    f_tmp.initialize();
    dvar_matrix F_fut_tmp(1,nfsh,1,nages);
    dvar_vector Ftottmp(1,nages);
		Ftottmp.initialize();
		F_fut_tmp.initialize();
		dvector p_lastyr(1,nfsh);

    dvar_matrix seltmp(1,nfsh,1,nages);
    seltmp.initialize();
    for (k=1;k<=nfsh;k++) seltmp(k) = (sel_fsh(k,endyr));

    for (k=1;k<=nfsh;k++) F_fut_tmp(k) = F(k,endyr);
    switch (iscenario)
    {
      case 1:
        // no multiplcation needed...it's 1.0...F_fut_tmp *= 1.087;
        // f_tmp = F35;
        // for (int k=1;k<=nfsh;k++) f_tmp(k) = SolveF2(endyr,nage_future(i), 1.0  * catch_lastyr(k));
        break;
      case 2:
        // f_tmp = SolveF2(endyr,nage_future(i), .75  * catch_lastyr );
        // for (int k=1;k<=nfsh;k++) f_tmp(k) = Fratio(k)*Fmsy; // mean(F(k,endyr));
        // for (int k=1;k<=nfsh;k++) f_tmp(k) = 0.75*mean(F(k,endyr));
        F_fut_tmp *= 0.75;
        break;
      case 3:
        // f_tmp = SolveF2(endyr,nage_future(i), 0.5 * catch_lastyr );
        //for (int k=1;k<=nfsh;k++) f_tmp(k) = .5*mean(F(k,endyr));
        // F_fut_tmp *= 1.25;
        for (int k=1;k<=nfsh;k++) F_fut_tmp(k) = seltmp(k)*Fratio(k)*Fmsy(s); // mean(F(k,endyr));
        break;
      case 4:
        // for (int k=1;k<=nfsh;k++) f_tmp(k) = .25*mean(F(k,endyr));
        F_fut_tmp *= 0.5;
        break;
      case 5:
        f_tmp = 0.0;
        F_fut_tmp = 0.0;
        break;
      /* case 6:
			// 15% increase...but doesn't seem to be working yet...
			  p_lastyr.initialize();
        for (int k=1;k<=nfsh;k++) p_lastyr(k) = catch_lastyr(k)/sum(catch_lastyr) ;
			  catch_lastyr = p_lastyr * 680. ; //591+0.15* 591.; OjO THIS NEEDS CHANGING AFTER 2019...
        f_tmp.initialize();
        f_tmp = SolveF2(endyr, catch_lastyr ,s);
        for (int k=1;k<=nfsh;k++) 
				  F_fut_tmp(k) = f_tmp(k)*seltmp(k);
				cout<<catch_lastyr<<endl;// <<F_fut_tmp<<endl;exit(1);
        break;
				*/
    }
    Z_future(s,i) = M(s,endyr);
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        // F_future(k,i) = sel_fsh(k,endyr) * F_fut_tmp(k);
        F_future(k,i)   = F_fut_tmp(k);
        Z_future(s,i)  += F_future(k,i);
      }
    }
    S_future(s,i) = mfexp(-Z_future(s,i));

FUNCTION Future_projections
  // Need to check on treatment of Fratio--whether it should be included or not
  SSB_fut.initialize();
  SSB_fut_1.initialize();
  SSB_fut_2.initialize();
  SSB_fut_3.initialize();
  SSB_fut_4.initialize();
  SSB_fut_5.initialize();
  // SSB_fut_6.initialize();
  catch_future.initialize();
  for (s=1;s<=nstk;s++)
  {
    for (int iscen=1;iscen<=5;iscen++)
    {
     // Future Sp_Biom set equal to estimated Sp_Biom w/ right lag
      // Sp_Biom_future(s)(styr_fut-rec_age,styr_fut-1) = Sp_Biom(s)(endyr-rec_age+1,endyr);
      for (i=styr_fut-rec_age;i<styr_fut;i++)
        Sp_Biom_future(s,i) = wt_mature(s) * elem_prod(natage(s,i),pow(S(s,i),spmo_frac)) ;

      nage_future(s,styr_fut)(2,nages) = ++elem_prod(natage(s,endyr)(1,nages-1),S(s,endyr)(1,nages-1));
      nage_future(s,styr_fut,nages)  += natage(s,endyr,nages)*S(s,endyr,nages);
      Sp_Biom_future(s,styr_fut)       = wt_mature(s) * elem_prod(nage_future(s,styr_fut),pow(S_future(s,styr_fut),spmo_frac)) ;
      // Future Recruitment (and Sp_Biom)
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(s,i,1)  = SRecruit( Sp_Biom_future(s,i-rec_age),cum_regs(s)+yy_sr(s,i) ) * mfexp(rec_dev_future(s,i)) ;     
        get_future_Fs(s,i,iscen);
        // Now graduate for the next year....
        nage_future(s,i+1)(2,nages) = ++elem_prod(nage_future(s,i)(1,nages-1),S_future(s,i)(1,nages-1));
        nage_future(s,i+1,nages)   += nage_future(s,i,nages)*S_future(s,i,nages);
        Sp_Biom_future(s,i) = wt_mature(s) * elem_prod(nage_future(s,i),pow(S_future(s,i),spmo_frac)) ;
      }
      nage_future(s,endyr_fut,1)  = SRecruit( Sp_Biom_future(s,endyr_fut-rec_age),cum_regs(s)+yy_sr(s,endyr_fut) ) * mfexp(rec_dev_future(s,endyr_fut)) ;     
      get_future_Fs(s,endyr_fut,iscen);
      Sp_Biom_future(s,endyr_fut)  = wt_mature(s) * elem_prod(nage_future(s,endyr_fut),pow(S_future(s,endyr_fut),spmo_frac)) ;
      if (iscen==1)
      {
        for (i=endyr+1;i<=endyr_fut;i++)
        {                   
          N_NoFsh(s,i,1)        = nage_future(s,i,1);
          // Adjustment for no-fishing recruits (ratio of R_nofish/R_fish)
          N_NoFsh(s,i,1)       *= SRecruit(Sp_Biom_NoFish(s,i-rec_age),cum_regs(s)+yy_sr(s,i)) / SRecruit(Sp_Biom_future(s,i-rec_age),cum_regs(s)+yy_sr(s,i));
          N_NoFsh(s,i)(2,nages) = ++N_NoFsh(s,i-1)(1,nages-1)*exp(-mean(natmort(s)));
          N_NoFsh(s,i,nages)   +=   N_NoFsh(s,i-1,nages)*exp(-mean(natmort(s)));
          Sp_Biom_NoFish(s,i)   = (N_NoFsh(s,i)*pow(exp(-mean(natmort(s))),spmo_frac) * wt_mature(s)); 
          // Sp_Biom_NoFishRatio(s,i)  = Sp_Biom_future(s,i) / Sp_Biom_NoFish(s,i) ;
        }
      }
      // Now get catch at future ages
      dvar_vector catage_tmp(1,nages);
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        catage_future(s,i).initialize();
        if (iscen!=5) 
        {
          for (k = 1 ; k<= nfsh ; k++)
          {
            if (sel_map(1,k) == s)
            {
              catage_tmp.initialize();
              catage_tmp = elem_prod(nage_future(s,i) , elem_prod(F_future(k,i) , 
                                    elem_div( ( 1.- S_future(s,i) ) , Z_future(s,i))));
              catage_future(s,i) += catage_tmp;
              catch_future(s,iscen,i)  += catage_tmp*wt_fsh(k,endyr);
            }
          }
        }
        SSB_fut(s,iscen,i) = Sp_Biom_future(s,i);
        switch (iscen)
        {
          case 1:
            SSB_fut_1(s,i) = Sp_Biom_future(s,i);
            break;
          case 2:
            SSB_fut_2(s,i) = Sp_Biom_future(s,i);
            break;
          case 3:
            SSB_fut_3(s,i) = Sp_Biom_future(s,i);
            break;
          case 4:
            SSB_fut_4(s,i) = Sp_Biom_future(s,i);
            break;
          case 5:
            SSB_fut_5(s,i) = Sp_Biom_future(s,i);
            break;
          /* case 6:
            SSB_fut_6(s,i) = Sp_Biom_future(s,i);
            break;
						*/
        }
      }
    }   //End of loop over F's
    Sp_Biom(s,endyr+1) = Sp_Biom_future(s,endyr+1);
  }

FUNCTION get_msy
 /*Function calculates used in calculating MSY and MSYL for a designated component of the
  population, given values for stock recruitment and selectivity...  
  Fmsy is the trial value of MSY example of the use of "funnel" to reduce the amount of storage for derivative calculations */

  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
    sumF(sel_map(1,k)) += sum(F(k,endyr));
  for (k=1;k<=nfsh;k++)
    Fratio(k) = sum(F(k,endyr)) / sumF(sel_map(1,k));

  dvar_vector Stmp(1,nstk);
  dvar_vector Rtmp(1,nstk);
  double df=1.e-05;
  for (s=1;s<=nstk;s++)
  {
    dvariable F1;
    F1.initialize();
    F1 = (0.8*natmortprior(mort_map(s,1))); //Ojo First year //nreg(s)
    dvariable F2;
    dvariable F3;
    dvariable yld1;
    dvariable yld2;
    dvariable yld3;
    dvariable dyld;
    dvariable dyldp;
    int breakout=0;
    // Newton Raphson stuff to go here
    for (int ii=1;ii<=8;ii++)
    {
      if (mceval_phase()&&(F1>5||F1<0.01)) 
      {
        ii=8;
        if (F1>5) F1=5.0; 
        else      F1=0.001; 
        breakout    = 1;
      }
      F2     = F1 + df*.5;
      F3     = F2 - df;
      // yld1   = yield(Fratio,F1, Stmp,Rtmp); // yld2   = yield(Fratio,F2,Stmp,Rtmp); // yld3   = yield(Fratio,F3,Stmp,Rtmp);
      yld1   = yield(Fratio,F1,s);
      yld2   = yield(Fratio,F2,s);
      yld3   = yield(Fratio,F3,s);
      dyld   = (yld2 - yld3)/df;                          // First derivative (to find the root of this)
      dyldp  = (yld2 + yld3 - 2.*yld1)/(.25*df*df);       // Second derivative (for Newton Raphson)
      if (breakout==0)
      {
        F1    -= dyld/dyldp;
      }
      else
      {
        if (F1>5) 
          cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
        else      
          cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
      }
    }
    {
      dvar_vector ttt(1,5);
      ttt         = yld(Fratio,F1,s);
      Fmsy(s)     = F1;
      Rtmp(s)     = ttt(3);
      MSY(s)      = ttt(2);
      Bmsy(s)     = ttt(1);
      MSYL(s)     = ttt(1)/Bzero(cum_regs(s)+yy_sr(s,endyr));
      lnFmsy(s)   = log(MSY(s)/ttt(5)); // Exploitation fraction relative to total biomass
      Bcur_Bmsy(s)= Sp_Biom(s,endyr)/Bmsy(s);

      dvariable FFtmp;
      FFtmp.initialize();
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == s)
          FFtmp += mean(F(k,endyr));
      Fcur_Fmsy(s)= FFtmp/Fmsy(s);
      Rmsy(s)     = Rtmp(s);
    }
  }

FUNCTION void get_msy(int iyr)
 /*Function calculates used in calculating MSY and MSYL for a designated component of the
  population, given values for stock recruitment and selectivity...  
  Fmsy is the trial value of MSY example of the use of "funnel" to reduce the amount of storage for derivative calculations */

  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
    sumF(sel_map(1,k)) += sum(F(k,iyr));
  for (k=1;k<=nfsh;k++)
    Fratio(k) = sum(F(k,iyr)) / sumF(sel_map(1,k));

  dvar_vector Stmp(1,nstk);
  dvar_vector Rtmp(1,nstk);
  double df=1.e-05;
  for (int istk=1;istk<=nstk;istk++)
  {
    dvariable F1;
    F1.initialize();
    F1 = (0.8*natmortprior(mort_map(istk,1))); //Ojo First year //nreg(s) //yy_sr(s,iyr)
    dvariable F2;
    dvariable F3;
    dvariable yld1;
    dvariable yld2;
    dvariable yld3;
    dvariable dyld;
    dvariable dyldp;
    int breakout=0;
    // Newton Raphson stuff to go here
    for (int ii=1;ii<=8;ii++)
    {
      if (mceval_phase()&&(F1>5||F1<0.01)) 
      {
        ii=8;
        if (F1>5) F1=5.0; 
        else      F1=0.001; 
        breakout    = 1;
      }
      F2     = F1 + df*.5;
      F3     = F2 - df;
      // yld1   = yield(Fratio,F1, Stmp,Rtmp); // yld2   = yield(Fratio,F2,Stmp,Rtmp); // yld3   = yield(Fratio,F3,Stmp,Rtmp);
      yld1   = yield(Fratio,F1,istk,iyr);
      yld2   = yield(Fratio,F2,istk,iyr);
      yld3   = yield(Fratio,F3,istk,iyr);
      dyld   = (yld2 - yld3)/df;                          // First derivative (to find the root of this)
      dyldp  = (yld2 + yld3 - 2.*yld1)/(.25*df*df);   // Second derivative (for Newton Raphson)
      if (breakout==0)
      {
        F1    -= dyld/dyldp;
      }
      else
      {
        if (F1>5) 
          cout<<"Fmsy v. high "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
        else      
          cout<<"Fmsy v. low "<< endl;// yld1<<" "<< yld2<<" "<< yld3<<" "<< F1<<" "<< F2<<" "<< F3<<" "<< endl;
      }
    }
    {
      dvar_vector ttt(1,5);
      ttt            = yld(Fratio,F1,istk,iyr);
      Fmsy(istk)     = F1;
      Rtmp(istk)     = ttt(3);
      MSY(istk)      = ttt(2);
      Bmsy(istk)     = ttt(1);
      MSYL(istk)     = ttt(1)/Bzero(cum_regs(istk)+yy_sr(istk,iyr));
      lnFmsy(istk)   = log(MSY(istk)/ttt(5)); // Exploitation fraction relative to total biomass
      Bcur_Bmsy(istk)= Sp_Biom(istk,iyr)/Bmsy(istk);

      dvariable FFtmp;
      FFtmp.initialize();
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == istk)
          FFtmp += mean(F(k,iyr));
      Fcur_Fmsy(istk)= FFtmp/Fmsy(istk);
      Rmsy(istk)     = Rtmp(istk);
    }
  }

FUNCTION dvar_vector yld(const dvar_vector& Fratio, const dvariable& Ftmp,int istk,int iyr)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector msy_stuff(1,5);
  dvariable phi;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  msy_stuff.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  seltmp.initialize();
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      seltmp(k) = sel_fsh(k,iyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = M(istk,iyr);
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);

  Ntmp(1) = 1.;
  for ( j=1 ; j < nages; j++ )
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/age class
  Ntmp(nages)  /= (1.- survtmp(nages)); 

  for (k=1;k<=nfsh;k++)
  {
    Ctmp.initialize();
    for ( j=1 ; j <= nages; j++ )
      Ctmp(j)      = Ntmp(j) * Fatmp(k,j) * (1. - survtmp(j)) / Ztmp(j);

    msy_stuff(2)  += wt_fsh(k,iyr) * Ctmp;
  }
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) ) * wt_mature(istk);
  // Req    = Requil(phi) * exp(sigmarsq/2);
  msy_stuff(5)  = Ntmp * wt_pop(istk);      
  msy_stuff(4)  = phi/phizero(cum_regs(istk)+yy_sr(istk,iyr)) ;       // SPR
  msy_stuff(3)  = Requil(phi,iyr,istk) ;       // Eq Recruitment
  msy_stuff(5) *= msy_stuff(3);       // BmsyTot
  msy_stuff(2) *= msy_stuff(3);       // MSY
  msy_stuff(1)  = phi*(msy_stuff(3)); // Bmsy
  RETURN_ARRAYS_DECREMENT();
  return msy_stuff;

 //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+ 
FUNCTION dvar_vector yld(const dvar_vector& Fratio, const dvariable& Ftmp, int istk)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector msy_stuff(1,5);
  dvariable phi;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  msy_stuff.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  seltmp.initialize();
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = M(istk,styr);
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);

  Ntmp(1) = 1.;
  for ( j=1 ; j < nages; j++ )
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/age class
  Ntmp(nages)  /= (1.- survtmp(nages)); 

  for (k=1;k<=nfsh;k++)
  {
    Ctmp.initialize();
    for ( j=1 ; j <= nages; j++ )
      Ctmp(j)      = Ntmp(j) * Fatmp(k,j) * (1. - survtmp(j)) / Ztmp(j);

    msy_stuff(2)  += wt_fsh(k,endyr) * Ctmp;
  }
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) ) * wt_mature(istk);
  // Req    = Requil(phi) * exp(sigmarsq/2);
  msy_stuff(5)  = Ntmp * wt_pop(istk);      
  msy_stuff(4)  = phi/phizero(cum_regs(istk)+yy_sr(istk,endyr)) ;       // SPR
  msy_stuff(3)  = Requil(phi,endyr,istk) ;       // Eq Recruitment
  msy_stuff(5) *= msy_stuff(3);       // BmsyTot
  msy_stuff(2) *= msy_stuff(3);       // MSY
  msy_stuff(1)  = phi*(msy_stuff(3)); // Bmsy
  RETURN_ARRAYS_DECREMENT();
  return msy_stuff;

FUNCTION dvariable yield(const dvar_vector& Fratio, const dvariable& Ftmp,int istk,int iyr)
  RETURN_ARRAYS_INCREMENT();
  dvariable phi;
  dvariable Req;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  dvariable   yield;
  yield.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      seltmp(k) = sel_fsh(k,iyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = M(istk,iyr);
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);

  Ntmp(1) = 1.;
  for ( j=1 ; j < nages; j++ )
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/age class
  Ntmp(nages)  /= (1.- survtmp(nages)); 

  for (k=1;k<=nfsh;k++)
  {
    Ctmp.initialize();
    for ( j=1 ; j <= nages; j++ )
      Ctmp(j)      = Ntmp(j) * Fatmp(k,j) * (1. - survtmp(j)) / Ztmp(j);

    yield  += wt_fsh(k,iyr) * Ctmp;
  }
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) )* wt_mature(istk);
  // Req    = Requil(phi) * mfexp(sigmarsq/2);
  Req    = Requil(phi,iyr,istk) ;
  yield *= Req;

  RETURN_ARRAYS_DECREMENT();
  return yield;

FUNCTION dvariable yield(const dvar_vector& Fratio, const dvariable& Ftmp,int istk)
  RETURN_ARRAYS_INCREMENT();
  dvariable phi;
  dvariable Req;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  dvariable   yield;
  yield.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  seltmp.initialize();
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = M(istk,styr);
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);

  Ntmp(1) = 1.;
  for ( j=1 ; j < nages; j++ )
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/age class
  Ntmp(nages)  /= (1.- survtmp(nages)); 

  for (k=1;k<=nfsh;k++)
  {
    Ctmp.initialize();
    for ( j=1 ; j <= nages; j++ )
      Ctmp(j)      = Ntmp(j) * Fatmp(k,j) * (1. - survtmp(j)) / Ztmp(j);

    yield  += wt_fsh(k,endyr) * Ctmp;
  }
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) )* wt_mature(istk);
  // Req    = Requil(phi) * mfexp(sigmarsq/2);
  Req    = Requil(phi,endyr,istk) ;
  yield *= Req;

  RETURN_ARRAYS_DECREMENT();
  return yield;

FUNCTION dvariable yield(const dvar_vector& Fratio, dvariable& Ftmp, dvariable& Stmp,dvariable& Req, int istk)
  RETURN_ARRAYS_INCREMENT();
  dvariable phi;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  dvariable   yield   = 0.;

  dvar_matrix seltmp(1,nfsh,1,nages);
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = M(istk,styr);
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);

  Ntmp(1) = 1.;
  for ( j=1 ; j < nages; j++ )
    Ntmp(j+1)  =   Ntmp(j) * survtmp(j); // Begin numbers in the next year/age class
  Ntmp(nages)  /= (1.- survtmp(nages)); 
  for (k=1;k<=nfsh;k++)
  {
    Ctmp.initialize();
    for ( j=1 ; j <= nages; j++ )
      Ctmp(j)      = Ntmp(j) * Fatmp(k,j) * (1. - survtmp(j)) / Ztmp(j);
    yield  += wt_fsh(k,endyr) * Ctmp;
  }
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) )* wt_mature(istk);
  // Req    = Requil(phi) * exp(sigmarsq/2);
  Req    = Requil(phi,endyr,istk) ;
  yield *= Req;
  Stmp   = phi*Req;

  RETURN_ARRAYS_DECREMENT();
  return yield;

FUNCTION Profile_F
  cout << "Doing a profile over F...."<<endl;
  ofstream prof_F("Fprof.yld");
 /* NOTE THis will need to be conditional on SrType too Function calculates 
  used in calculating MSY and MSYL for a designated component of the
  population, given values for stock recruitment and selectivity...  
  Fmsy is the trial value of MSY example of the use of "funnel" to 
  reduce the amount of storage for derivative calculations */
  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
    sumF(sel_map(1,k)) += sum(F(k,endyr));
  for (k=1;k<=nfsh;k++)
    Fratio(k) = sum(F(k,endyr)) / sumF(sel_map(1,k));
  dvariable Stmp;
  dvariable Rtmp;
  double df=1.e-7;
  dvariable F1=.05;
  dvariable F2;
  dvariable F3;
  dvariable yld1;
  dvariable yld2;
  dvariable yld3;
  dvariable dyld;
  dvariable dyldp;
  prof_F <<"Profile of stock, yield, and recruitment over F"<<endl;
  prof_F << model_name<<" "<<datafile_name<<endl;
  prof_F <<endl<<endl<<"F  Stock  Yld  Recruit SPR"<<endl;
  for (s=1;s<=nstk;s++)
  {
    for (r=1;r<=nreg(s);r++)
      prof_F <<0.0<<" "<< Bzero(cum_regs(s)+r) <<" "<<0.0<<" "<<Rzero(cum_regs(s)+r)<< " 1.00"<<endl; //Ojo
    dvar_vector ttt(1,5);
    for (int ii=1;ii<=500;ii++)
    {
      F1    = double(ii)/500;
      yld1  = yield(Fratio,F1,Stmp,Rtmp,s);
      ttt   = yld(Fratio,F1,s);
      prof_F <<F1<<" "<< ttt << endl; 
    }
    prof_F <<endl;
  }

FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp, const int& Nsr_tmp)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
  switch (SrType)
  {
    case 1:
      RecTmp = elem_prod((Stmp / phizero(Nsr_tmp)) , mfexp( alpha(Nsr_tmp) * ( 1. - Stmp / Bzero(Nsr_tmp) ))) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = elem_prod(Stmp , 1. / ( alpha(Nsr_tmp) + beta(Nsr_tmp) * Stmp));        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec(Nsr_tmp));                    //Avg recruitment
      break;
    case 4:
      RecTmp = elem_prod(Stmp , mfexp( alpha(Nsr_tmp)  - Stmp * beta(Nsr_tmp))) ; //Old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(const double& Stmp, const int& Nsr_tmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero(Nsr_tmp)) * mfexp( alpha(Nsr_tmp) * ( 1. - Stmp / Bzero(Nsr_tmp) )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha(Nsr_tmp) + beta(Nsr_tmp) * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec(Nsr_tmp));                    //Avg recruitment
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha(Nsr_tmp)  - Stmp * beta(Nsr_tmp)) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(const dvariable& Stmp,const int& Nsr_tmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero(Nsr_tmp)) * mfexp( alpha(Nsr_tmp) * ( 1. - Stmp / Bzero(Nsr_tmp) )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha(Nsr_tmp) + beta(Nsr_tmp) * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec(Nsr_tmp) );                    //Avg recruitment
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha(Nsr_tmp)  - Stmp * beta(Nsr_tmp)) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

 //=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
FUNCTION Get_Bzero
  Bzero.initialize();
  Rzero    =  mfexp(log_Rzero); 

  dvar_matrix survtmp(1,nstk,1,nages);
  for (s=1;s<=nstk;s++)
    survtmp(s) = mfexp(-M(s,styr));

  dvar3_array natagetmp(1,nstk,styr_rec,styr,1,nages);
  natagetmp.initialize();

  Sp_Biom.initialize();
  for (s=1;s<=nstk;s++)
  {
    natagetmp(s,styr_rec,1) = Rzero(cum_regs(s)+1);
    for (j=2; j<=nages; j++)
      natagetmp(s,styr_rec,j) = natagetmp(s,styr_rec,j-1) * survtmp(s,j-1);
    natagetmp(s,styr_rec,nages) /= (1.-survtmp(s,nages));
    
    Bzero(cum_regs(s)+1) = elem_prod(wt_mature(s) , pow(survtmp(s),spmo_frac))*natagetmp(s,styr_rec);
    Sp_Biom(s)(styr_sp,styr_rec-1) = Bzero(cum_regs(s)+1);
    for (i=styr_rec;i<styr;i++)
    {
      Sp_Biom(s,i) = elem_prod(natagetmp(s,i),pow(survtmp(s),spmo_frac)) * wt_mature(s);
      // natagetmp(s,i,1)          = mfexp(rec_dev(s,i) + log_Rzero(cum_regs(s)+1)); // OjO numbers a function of mean not SR curve...
      natagetmp(s,i,1)          = mfexp(rec_dev(s,i) + mean_log_rec(cum_regs(s)+yy_sr(s,i)));
      natagetmp(s,i+1)(2,nages) = ++elem_prod(natagetmp(s,i)(1,nages-1),mfexp(-M(s,styr)(1,nages-1)) );
      natagetmp(s,i+1,nages)   += natagetmp(s,i,nages)*mfexp(-M(s,styr,nages));
    }
    // This sets first year recruitment as deviation from mean recruitment (since SR curve can
    // be defined for different periods and is treated semi-independently)
    natagetmp(s,styr,1)   = mfexp(rec_dev(s,styr) + mean_log_rec(cum_regs(s)+yy_sr(s,styr)));
    mod_rec(s)(styr_rec,styr) = column(natagetmp(s),1);
    natage(s,styr)  = natagetmp(s,styr); // OjO
    Sp_Biom(s,styr) = elem_prod(natagetmp(s,styr),pow(survtmp(s),spmo_frac)) * wt_mature(s); 
  }

FUNCTION dvariable Requil(dvariable& phi, int iyr, int istk)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  int ireg=cum_regs(istk)+yy_sr(istk,iyr);
  switch (SrType)
  {
    case 1:
      RecTmp =  Bzero(ireg) * (alpha(ireg) + log(phi) - log(phizero(ireg)) ) / (alpha(ireg)*phi);
      break;
    case 2:
      RecTmp =  (phi-alpha(ireg))/(beta(ireg)*phi);
      break;
    case 3:
      RecTmp =  mfexp(mean_log_rec(ireg));
      break;
    case 4:
      RecTmp =  (log(phi)+alpha(ireg)) / (beta(ireg)*phi); //RecTmp =  (log(phi)/alpha + 1.)*beta/phi;
      break;
  }
  // Req    = Requil(phi) * exp(sigmarsq/2);
  // return RecTmp* exp(sigmarsq/2);
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;
  
FUNCTION write_mceval_hdr
    for (k=1;k<=nind;k++)
      mceval<< " q_ind_"<< k<< " ";
    mceval<<"M steepness depletion MSY MSYL Fmsy Fcur_Fmsy Bcur_Bmsy Bmsy totbiom_"<<endyr<<" "<< 
    " F35          "<< 
    " F40          "<< 
    " F50          "<< 
    " fut_SPB_Fmsy_"<< endyr_fut<<" "<< 
    " fut_SPB_F50%_"<< endyr_fut<<" "<< 
    " fut_SPB_F40%_"<< endyr_fut<<" "<< 
    " fut_SPB_F35%_"<< endyr_fut<<" "<< 
    " fut_SPB_F0_"  << endyr_fut<<" "<< 
    " fut_catch_Fmsy_"<<styr_fut<<" "<<  
    " fut_catch_F50%_"<<styr_fut<<" "<<  
    " fut_catch_F40%_"<<styr_fut<<" "<<  
    " fut_catch_F35%_"<<styr_fut<<" "<<  endl;

//+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+ 
REPORT_SECTION
  if (last_phase())
  {
    save_gradients(gradients);
    int nvar1=initial_params::nvarcalc(); // get the number of active parameters
    int ndvar=stddev_params::num_stddev_calc();
    int offset=1;
    dvector param_values(1,nvar1+ndvar);
    initial_params::copy_all_values(param_values,offset);
    for (int i=0;i<initial_params::num_initial_params;i++)
    {
      if (withinbound(0,(initial_params::varsptr[i])->phase_start, initial_params::current_phase))
      {
        int sc = (initial_params::varsptr[i])->size_count();
        if (sc>0)
        {
          // write_input_log << "# " << initial_params::varsptr[i]->label() << endl<<param_values(i)<<"\n" << endl; 
        } 
      }
    }
        
    for (k=1;k<=nfsh;k++)
      if (nyrs_fsh_age(k)>0)
        Francis_wts(k) = calc_Francis_weights(oac_fsh(k), eac_fsh(k), n_sample_fsh_age(k) );

      // age_like_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(eac_ind(k,i) + 0.001 ) ;
    for (k=1;k<=nind;k++)
    {
      // age_like_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(eac_ind(k,i) + 0.001 ) ;
      if (nyrs_ind_age(k)>0)
        Francis_wts(k+nfsh) = calc_Francis_weights(oac_ind(k), eac_ind(k), n_sample_ind_age(k) );
    }

    if (!Popes)
      for (k=1;k<=nfsh;k++)
        Ftot(sel_map(1,k)) += F(k);
    for (int r=1;r<=nmort;r++)
    {
      log_param(Mest(r));
    }
    //log_param(mean_log_rec);
    for (int r=1;r<=nregs;r++)
    {
      log_param(mean_log_rec(r));
    }
    for (int r=1;r<=nrec;r++)
    {
      log_param(steepness(r));
    }
    for (int r=1;r<=nregs;r++)
    {
      log_param(log_Rzero(r));
    }
    //log_param(log_Rzero);
    log_param(rec_dev);
    for (int r=1;r<=nrec;r++)
    {
      log_param(log_sigmar(r));
    }
    log_param(fmort);
    // log_param(log_selcoffs_fsh);
    // log_param(log_sel_spl_fsh);
    // log_param(logsel_slope_fsh);
    // log_param(sel50_fsh);
    // log_param(logsel_dslope_fsh);
    // log_param(seld50_fsh);
    log_param(rec_dev_future);
    // log_param(log_q_ind);
    // log_param(log_q_power_ind);
    // log_param(log_selcoffs_ind);
    // log_param(logsel_slope_ind);
    // log_param(logsel_dslope_ind);
    // log_param(sel50_ind);
    // log_param(seld50_ind);
  }
    
  if (oper_mod)
    Oper_Model();

  cout <<"==============================================================="<<endl;
  if(last_phase())
    cout<<"||  ++++++ Completed phase "<<current_phase()<<" In last phase now +++++"<< endl<<"||"<<endl<<"||  "<<cntrlfile_name <<endl;
  else
    cout<<"||  ++++++ Completed phase "<<current_phase()<<" ++++++++++++++++"<< endl<<"||"<<endl<<"||  "<<cntrlfile_name <<endl;
  cout<<"||"<<endl<<"||"<<endl;
  cout <<"_______________________________________________________________"<<endl;
    adstring comma = adstring(","); 
    report << model_name<<" "<< endl<< endl;
    report << "Estimated annual F's " << endl;
    Fmort.initialize();
    for (k=1;k<=nfsh;k++)
      for (i=styr;i<=endyr;i++) 
        Fmort(sel_map(1,k),i) += mean(F(k,i));
    report << Fmort<<endl;
    report << "Total mortality (Z)"<<endl;
    for (s=1;s<=nstk;s++)
      report << Z(s)<<endl;
    report << "Estimated numbers of fish " << endl;
    for (s=1;s<=nstk;s++)
      for (i=styr;i<=endyr;i++) 
        report <<"       Year: "<< i << " "<< natage(s,i) << endl;
    report << endl<< "Estimated F mortality " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "Fishery "<< k <<" : "<< endl ;
      for (i=styr;i<=endyr;i++) 
        report << "        Year: "<<i<<" "<<F(k,i)<<  " "<< endl;
    }

    report << endl<< "survey q " << endl;
    report <<q_ind<<endl;
    report << endl<< "Observed survey values " << endl;
    for (k=1;k<=nind;k++)
    {
      int ii=1;
      int istk=sel_map(1,k+nfsh);
      report <<endl<< "Yr_Obs_Pred_Survey "<< k <<" : "<< endl ;
      for (int iyr=styr;iyr<=endyr;iyr++)
      {
        dvariable pred_tmp ;
        if (ii<=nyrs_ind(k))
        {
          pred_tmp = q_ind(k,ii) * pow(elem_prod(natage(istk,iyr),pow(S(istk,iyr),ind_month_frac(k))) * 
                        elem_prod(sel_ind(k,iyr) , wt_ind(k,iyr)),q_power_ind(k));
          if (yrs_ind(k,ii)==iyr)
          {
            report << iyr<< " "<< 
                     obs_ind(k,ii) << " "<< pred_tmp <<endl;
            ii++;
          }
          else
            report << iyr<< " -1 "<< " "<< pred_tmp   <<endl;
        }
      }
    }

    report << endl<< "Survey_Q:  "<<q_ind << endl;

    report << endl<< "Observed Prop " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "ObsFishery "<< k <<" : "<< endl ;
      for (i=1;i<=nyrs_fsh_age(k);i++) 
        report << yrs_fsh_age(k,i)<< " "<< oac_fsh(k,i) << endl;
    }
    for (k=1;k<=nfsh;k++)
    {
      report << "Pobs_length_fishery_"<< (k) <<""<< endl;
      for (i=1;i<=nyrs_fsh_length(k);i++) 
        report << yrs_fsh_length(k,i)<< " "<< olc_fsh(k,i) << endl;
      report   << endl;
    }
    report << endl<< "Predicted prop  " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "PredFishery "<< k <<" : "<< endl;
      for (i=1;i<=nyrs_fsh_age(k);i++) 
        report << yrs_fsh_age(k,i)<< " "<< eac_fsh(k,i) << endl;
    }
    for (k=1;k<=nfsh;k++)
    {
      report << "Pred_length_fishery_"<< (k) <<""<< endl;
      for (i=1;i<=nyrs_fsh_length(k);i++) 
        report << yrs_fsh_length(k,i)<< " "<< elc_fsh(k,i) << endl;
      report   << endl;
    }
    report << endl<< "Observed prop Survey" << endl;
    for (k=1;k<=nind;k++)
    {
      report << "ObsSurvey "<<k<<" : "<<  endl;
      for (i=1;i<=nyrs_ind_age(k);i++) 
        report << yrs_ind_age(k,i)<< " "<< oac_ind(k,i) << endl;
    }
    for (k=1;k<=nind;k++)
    {
      report << "Pobs_length_survey "<< (k) <<""<<  endl;
      for (i=1;i<=nyrs_ind_length(k);i++) 
        report << yrs_ind_length(k,i)<< " "<< olc_ind(k,i) << endl;
    }
    report << endl<< "Predicted prop Survey" << endl;
    for (k=1;k<=nind;k++)
    {
      report << "PredSurvey "<<k<<" : "<<  endl;
      for (i=1;i<=nyrs_ind_age(k);i++) 
        report << yrs_ind_age(k,i)<< " "<< eac_ind(k,i) << endl;
    }
    for (k=1;k<=nind;k++)
    {
      report << "Pred_length_survey "<< (k) <<""<<  endl;
      for (i=1;i<=nyrs_ind_length(k);i++) 
        report << yrs_ind_length(k,i)<< " "<< elc_ind(k,i) << endl;
    }
    report << endl<< "Observed catch biomass " << endl;
    report << catch_bio << endl;
    report << "predicted catch biomass " << endl;
    report << pred_catch << endl;

    report << endl<< "Estimated annual fishing mortality " << endl;
    for (k=1;k<=nfsh;k++)
      report << " Average_F_Fshry_"<<k<< " Full_selection_F_Fshry_"<<k;

    report << endl;
    for (i=styr;i<=endyr;i++)
    {
      report<< i<< " ";
      for (k=1;k<=nfsh;k++)
        report<< mean(F(k,i)) <<" "<< mean(F(k,i))*max(sel_fsh(k,i)) << " ";

      report<< endl;
    }
    report << endl<< "Selectivity" << endl;
    for (k=1;k<=nfsh;k++)
      for (i=styr;i<=endyr;i++)
        report << "Fishery "<< k <<"  "<< i<<" "<<sel_fsh(k,i) << endl;
    for (k=1;k<=nind;k++)
      for (i=styr;i<=endyr;i++)
        report << "Survey  "<< k <<"  "<< i<<" "<<sel_ind(k,i) << endl;

    report << endl<< "Stock Recruitment stuff "<< endl;
    for (s=1;s<=nstk;s++)
      for (i=styr_rec;i<=endyr;i++)
        if (active(log_Rzero(cum_regs(s)+yy_sr(s,i))))
          report << i<< " "<<Sp_Biom(s,i-rec_age)<< " "<< SRecruit(Sp_Biom(s,i-rec_age),cum_regs(s)+yy_sr(s,i))<< " "<< mod_rec(s,i)<<endl;
        else 
          report << i<< " "<<Sp_Biom(s,i-rec_age)<< " "<< " 999" << " "<< mod_rec(s,i)<<endl;

    report << endl<< "Curve to plot "<< endl;
    report <<"stock Recruitment"<<endl;
    for (r=1;r<=nregs;r++)
    {
      report <<"0 0 "<<endl;
      dvariable stock;
      for (i=1;i<=300;i++)
      {
        stock = double (i) * Bzero(r) /250.;
        if (active(log_Rzero(r)))
          report << stock <<" "<< SRecruit(stock, r)<<endl;
        else
          report << stock <<" 99 "<<endl;
      }
    }

    report   << endl<<"Likelihood Components" <<endl;
    report   << "----------------------------------------- " <<endl;
    report   << "  catch_like age_like_fsh length_like_fsh sel_like_fsh ind_like age_like_ind length_like_ind sel_like_ind rec_like fpen post_priors_indq post_priors residual total"<<endl;
    report   << " "<<obj_comps<<endl;

    obj_comps(13)= obj_fun - sum(obj_comps(1,12)) ; // Residual 
    obj_comps(14)= obj_fun ;                  // Total
    report   <<"  catch_like       "<<setw(10)<<obj_comps(1) <<endl
             <<"  age_like_fsh     "<<setw(10)<<obj_comps(2) <<endl
             <<"  length_like_fsh  "<<setw(10)<<obj_comps(3) <<endl
             <<"  sel_like_fsh     "<<setw(10)<<obj_comps(4) <<endl
             <<"  ind_like         "<<setw(10)<<obj_comps(5) <<endl
             <<"  age_like_ind     "<<setw(10)<<obj_comps(6) <<endl
             <<"  length_like_ind  "<<setw(10)<<obj_comps(7) <<endl
             <<"  sel_like_ind     "<<setw(10)<<obj_comps(8) <<endl
             <<"  rec_like         "<<setw(10)<<obj_comps(9) <<endl
             <<"  fpen             "<<setw(10)<<obj_comps(10)<<endl
             <<"  post_priors_indq "<<setw(10)<<obj_comps(11)<<endl
             <<"  post_priors      "<<setw(10)<<obj_comps(12)<<endl
             <<"  residual         "<<setw(10)<<obj_comps(13)<<endl
             <<"  total            "<<setw(10)<<obj_comps(14)<<endl;
    report   << endl;
    report   << "Fit to Catch Biomass "<<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nfsh;k++)
      report << "  Catch_like_Fshry_#"<< k <<"  "<< catch_like(k) <<endl;
    report   << endl;

    report << "Age likelihoods for fisheries :"<<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nfsh;k++)
      report << "  Age_like_Fshry_#"<< k <<"  "<< age_like_fsh(k) <<endl;
    report   << endl;

    report   << "Selectivity penalties for fisheries :"<<endl;
    report   << "-------------------------" <<endl;
    report   << "  Fishery Curvature_Age Change_Time Dome_Shaped"<<endl;
    for (k=1;k<=nfsh;k++)
      report << "  Sel_Fshry_#"<< k <<"  "<< sel_like_fsh(k,1) <<" "<<sel_like_fsh(k,2)<<" "<<sel_like_fsh(k,3)<< endl;
    report   << endl;
  
    report   << "survey Likelihood(s) " <<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nind;k++)
      report << "  Survey_Index_#"<< k <<"  " << ind_like(k)<<endl;
    report   << endl;

    report << setw(10)<< setfixed() << setprecision(5) <<endl;
    report   << "Age likelihoods for surveys :"<<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nind;k++)
      report << "  Age_Survey_#"<< k <<"  " << age_like_ind(k)<<endl;
    report   << endl;

    report   << "Selectivity penalties for surveys :"<<endl;
    report   << "-------------------------" <<endl;
    report   << "  Survey Curvature_Age Change_Time Dome_Shaped"<<endl;
    for (k=1;k<=nind;k++)
      report << "  Sel_Survey_#"<< k <<"  "<< sel_like_ind(k,1) <<" "<<sel_like_ind(k,2)<<" "<<sel_like_ind(k,3)<< endl;
    report   << endl;

    report << setw(10)<< setfixed() << setprecision(5) <<endl;
    report   << "Recruitment penalties: " <<endl<<rec_like<<endl;
    report   << "-------------------------" <<endl;
    report   << "  (sigmar)            " <<sigmar<<endl;
    report   << "  S-R_Curve           " <<column(rec_like,1)<< endl;
    report   << "  Regularity          " <<column(rec_like,2)<< endl;
    report   << "  Future_Recruits     " <<column(rec_like,3)<< endl;
    report   << endl;

    report   << "F penalties:          " <<endl;
    report   << "-------------------------" <<endl;
    report   << "  Avg_F               " <<fpen(1) <<endl;
    report   << "  Effort_Variability  " <<fpen(2) <<endl;
    report   << endl;

    report   << "Contribution of Priors:"<<endl;
    report   << "-------------------------" <<endl;
    report   << "Source                ";
    report   <<           " Posterior";
    report   <<           " Param_Val";
    report   <<           " Prior_Val";
    report   <<           "  CV_Prior"<<endl;
  // (*ad_printf)("f = %lf\n",value(f));
    for (k=1;k<=nind;k++)
    {
      report << "Q_Survey_#"<< k <<"           "
             << setw(10)<<post_priors_indq(k) 
             << setw(10)<< q_ind(k)
             << setw(10)<< qprior(k)
             << setw(10)<< cvqprior(k)<<endl;

      report << "Q_power_Survey_#"<< k <<"           "
             << setw(10)<<post_priors_indq(k) 
             << setw(10)<< q_power_ind(k)
             << setw(10)<< q_power_prior(k)
             << setw(10)<< cvq_power_prior(k)<<endl;
    }

    // writerep(post_priors(1),repstring);
    report   << "Natural_Mortality     "
             << setw(10)<< column(post_priors,1)(1,nmort)
             << setw(10)<< M
             << setw(10)<< natmortprior
             << setw(10)<< cvnatmortprior <<endl;
    report   << "Steepness             "
             << setw(10)<< column(post_priors,2)(1,nrec)
             << setw(10)<< steepness
             << setw(10)<< steepnessprior
             << setw(10)<< cvsteepnessprior <<endl;
    report   << "SigmaR                "
             << setw(10)<< column(post_priors,3)(1,nrec)
             << setw(10)<< sigmar
             << setw(10)<< sigmarprior
             << setw(10)<< cvsigmarprior <<endl;
    report   << endl;
    report<<"Num_parameters_Estimated "<<initial_params::nvarcalc()<<endl;
    
  report <<cntrlfile_name<<endl;
  report <<datafile_name<<endl;
  report <<model_name<<endl;
  if (SrType==2) 
    report<< "Beverton-Holt" <<endl;
  else
    report<< "Ricker" <<endl;
  report<<"Steepnessprior,_CV,_phase: " <<steepnessprior<<" "<<
    cvsteepnessprior<<" "<<
    phase_srec<<" "<< endl;

  report<<"sigmarprior,_CV,_phase: " <<sigmarprior<<" "<<  cvsigmarprior <<" "<<phase_sigmar<<endl;

  report<<"Rec_estimated_in_styr_endyr: " <<styr_rec    <<" "<<endyr        <<" "<<endl;
  for (r=1;r<=nregs;r++)
  {
    int istk = stk_reg_map(1,r);
    int ireg = stk_reg_map(2,r);
    report<<"SR_Curve_fit__in_styr_endyr_" <<r<<" : " <<styr_rec_est(istk,ireg)<<" "<<endyr_rec_est(istk,ireg)<<" "<<endl;
  }
  report<<"Model_styr_endyr:            " <<styr        <<" "<<endyr        <<" "<<endl;

  report<<"M_prior,_CV,_phase "<< natmortprior<< " "<< cvnatmortprior<<" "<<phase_M<<endl;
  report<<"qprior,_CV,_phase " <<qprior<<" "<<cvqprior<<" "<< phase_q<<endl;
  report<<"q_power_prior,_CV,_phase " <<q_power_prior<<" "<<cvq_power_prior<<" "<< phase_q_power<<endl;

  report<<"cv_catchbiomass: " <<cv_catchbiomass<<" "<<endl;
  report<<"Projection_years "<< nproj_yrs<<endl;
  for (k=1;k<=nfsh;k++)
    report << "Fsh_sel_opt_fish: "<<k<<" "<<fsh_sel_opt(k)<<" "<<sel_change_in_fsh(k)<<endl;
  for (k=1;k<=nind;k++)
    report<<"Survey_Sel_Opt_Survey: " <<k<<" "<<(ind_sel_opt(k))<<endl;
    
  report <<"Phase_survey_Sel_Coffs: "<<phase_selcoff_ind<<endl; 
  report <<"Fshry_Selages: " << nselages_in_fsh  <<endl;
  report <<"Survy_Selages: " << nselages_in_ind <<endl;
  report << "Phase_for_age-spec_fishery "<<phase_selcoff_fsh<<endl;
  report << "Phase_for_logistic_fishery "<<phase_logist_fsh<<endl;
  report << "Phase_for_dble_logistic_fishery "<<phase_dlogist_fsh<<endl;
  report << "Phase_for_age-spec_survey  "<<phase_selcoff_ind<<endl;
  report << "Phase_for_logistic_survey  "<<phase_logist_ind<<endl;
  report << "Phase_for_dble_logistic_indy "<<phase_dlogist_ind<<endl;

  for (k=1; k<=nfsh;k++)
  {
    report <<"Number_of_select_changes_fishery: "<<k<<" "<<n_sel_ch_fsh(k)<<endl;
    report<<"Yrs_fsh_sel_change: "<<yrs_sel_ch_fsh(k)<<endl;
    report << "sel_change_in: "<<sel_change_in_fsh(k) << endl;
  }
  for (k=1; k<=nind;k++)
  {
    report <<"Number_of_select_changes_survey: "<<k<<" "<<n_sel_ch_ind(k)<<endl;
    report<<"Yrs_ind_sel_change: "<<yrs_sel_ch_ind(k)<<endl;
    report << "sel_change_in: "<<sel_change_in_ind(k) << endl;
  }

FUNCTION write_msy_out
  ofstream msyout("msyout.dat");
  msyout << " # Natural Mortality       " <<endl;
  for (s=1;s<=nstk;s++) 
    msyout <<M(s) <<" ";
  msyout <<endl;
  msyout << spawnmo<< "  # Spawnmo                   " <<endl;
  msyout <<"# Wt spawn"<<endl<< wt_pop<< endl;
  msyout <<"# Wt fish"<<endl;
  for (k=1;k<=nfsh;k++) 
    msyout <<wt_fsh(k,endyr)<< " ";
  msyout <<endl;
  msyout <<"# Maturity"<<endl<< maturity<< endl;
  msyout <<"# selectivity"<<endl;
  for (k=1;k<=nfsh;k++) 
    msyout<< sel_fsh(k,endyr) <<" ";
  msyout<< endl;
  msyout<<"Srec_Option "<<SrType<< endl;
  msyout<<"Alpha "<<alpha<< endl;
  msyout<<"beta "<<beta<< endl;
  msyout<<"steepness "<<steepness<< endl;
  msyout<<"Bzero "<<Bzero<< endl;
  msyout<<"Rzero "<<Rzero<< endl;

FUNCTION write_projout
// Function to write out data file for projection model....
  ofstream projout( projfile_name );
  
  projout <<"# "<<model_name <<" "<< projfile_name<<endl;
  projout <<"123  # seed"<<endl;
  // Flag to tell if this is a SSL species                 
  projout << "1 # Flag to tell if this is a SSL forage species                 "<<endl;
  projout << "0 # Flag to Dorn's version of a constant buffer                  "<<endl;
  // Flag to solve for F in first year or not 0==don't solve
  projout<< " 1 # Flag to solve for F in first year or not 0==don't solve"<<endl;
  // Flag to use 2nd-year catch/TAC
  projout<< "0 # Flag to use 2nd-year catch/TAC"<<endl;
  projout << nstk<<"   # Number of stocks"<<endl;
  projout << nfsh<<"   # Number of fisheries"<<endl;
  projout <<"14   # Number of projection years"<<endl;
  projout <<"1000 # Number of simulations"<<endl;
  projout <<endyr<< " # Begin year of projection" <<endl;
  projout <<nages<< " # Number of ages" <<endl;
  for (s=1;s<=nstk;s++) 
    projout <<M(s) <<" "<<endl;
  projout << " # Natural Mortality       " <<endl;
  double sumtmp;
  sumtmp = 0.;
  for (k=1;k<=nfsh;k++) 
    sumtmp += catch_bio(k,endyr);
  projout << sumtmp<< " # TAC in current year (assumed catch) " <<endl;
  projout << sumtmp<< " # TAC in current year+1 (assumed catch) " <<endl;
  for (k=1;k<=nfsh;k++) 
    projout <<  F(k,endyr)/mean((F(k,endyr)))<<" "<<endl;
   //  + fmort_dev(k,endyr)) /Fmort(endyr)<<" ";

  projout << "   # Fratio                  " <<endl;
  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
  {
    Fratio(k) = sum(F(k,endyr)) ;
    sumF(sel_map(1,k)) += Fratio(k) ;
  }
  for (k=1;k<=nfsh;k++)
    Fratio(k) /= sumF(sel_map(1,k));
  projout << Fratio         <<endl;
  projout <<"  # average f" <<endl;
  projout << " 1  # author f                  " <<endl;
  projout << spawnmo<< "  # Spawnmo                   " <<endl;
  projout <<"# Wt spawn"<<endl<< wt_pop<< endl;
  projout <<"# Wt fish"<<endl;
  for (k=1;k<=nfsh;k++) 
    projout <<wt_fsh(k,endyr)<< " ";
  projout <<endl;
  projout <<"# Maturity"<<endl<< maturity<< endl;
  projout <<"# selectivity"<<endl;
  for (k=1;k<=nfsh;k++) 
    projout<< sel_fsh(k,endyr) <<" "<<endl;
  projout<< endl;
  projout <<"# natage"<<endl;
  for (s=1;s<=nstk;s++)
    projout << natage(s,endyr) << endl;
  if (styr<(1977-rec_age-1))
  {
    projout <<"#_N_recruitment_years (not including last 1 estimates)"<<endl<<endyr-(1977+rec_age+1) << endl;
    projout <<"#_Recruitment_start_at_1977_yearclass=1978_for_age_1_recruits"<<yy(1977+rec_age,endyr-1)<<endl;
    for (s=1;s<=nstk;s++)
      projout << mod_rec(s)(1977+rec_age,endyr-1)<< endl;
  }

FUNCTION write_proj
 ofstream newproj("proj.dat");
// Function to write out data file for new Ianelli 2005 projection model....
 newproj <<"#Species name here:"<<endl;
 newproj <<model_name+"_"+datafile_name<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 // Need to correct for maxf standardization 

 dvector seltmp(1,nages);
 dvar_vector sumF(1,nstk);
 sumF.initialize();
 seltmp.initialize();
 for (k=1;k<=nfsh;k++)
 {
   Fratio(k) = sum(F(k,endyr)) ;
   sumF(sel_map(1,k)) += value(Fratio(k)) ;
 }
 for (k=1;k<=nfsh;k++)
   Fratio(k) /= sumF(sel_map(1,k));
 // compute a 5-year recent average fishery-aggregated selectivity for output to projection model
 for (k=1;k<=nfsh;k++)
   for (j=1;j<=nages;j++)
     seltmp(j) += value(Fratio(k))*(value(sel_fsh(k,endyr,j)) 
                 +value(sel_fsh(k,endyr-1,j))  
                 +value(sel_fsh(k,endyr-2,j))  
                 +value(sel_fsh(k,endyr-3,j))  
                 +value(sel_fsh(k,endyr-4,j))
                 )/5.;  

 for (s=1;s<=nstk;s++)
   newproj << mean(Fmort(s)(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << spmo_frac*12+1<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#_Natural_Mortality" << aa << endl;
   for (s=1;s<=nstk;s++) newproj <<M(s,endyr)<<" "; newproj<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<2.*maturity<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt_pop<< endl;
 newproj <<"#_Wt_at_age_fishery" <<aa<<endl;
   for (k=1;k<=nfsh;k++) newproj<<wt_fsh(k,endyr) << endl;
 newproj <<"#" <<endl;

 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl;
 seltmp = value(sel_fsh(1,endyr)) +value(sel_fsh(1,endyr-1))  +value(sel_fsh(1,endyr-2));  
 newproj << seltmp/max(seltmp)<<endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl;
   for (s=1;s<=nstk;s++) newproj <<natage(s,endyr)<< endl;
  if (styr<=1977)
  {
   newproj <<"#_N_recruitment_years (not including last estimate)"<<endl<<endyr-(1977+rec_age) << endl;
   newproj <<"#_Recruitment_start_at_1977_yearclass=1978_for_age_1_recruits"<<yy(1977+rec_age,endyr-1)<<endl;
     for (s=1;s<=nstk;s++) newproj <<mod_rec(s)(1977+rec_age,endyr-1)<< endl;
  }

 newproj <<"#_Spawning biomass "<<endl;
   for (s=1;s<=nstk;s++) newproj <<Sp_Biom(s)(styr-rec_age,endyr-rec_age)/1000<< endl;
 newproj.close();
 
RUNTIME_SECTION
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-5,1e-5
  maximum_function_evaluations 100,100,200,300,2500

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(4500);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(4500);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize=500000000;

FINAL_SECTION
  // Calc_Dependent_Vars();
  write_proj();
  write_projout();
  // write_msy_out();
  Profile_F();
  Write_R();
  
GLOBALS_SECTION
  #include <admodel.h>
	#undef truth
	#define truth(object) trudat << #object "\n" << object << endl;
	#undef REPORT
	#define REPORT(object) REPORT << #object "\n" << object << endl;
	#undef R_Report
	#define R_Report(object) R_report << "$"#object "\n" << object << endl;
	/**
	 \def log_input(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef log_input
	#define log_input(object) write_input_log << "# " #object "\n" << object << endl;
	#undef log_param
  // #define log_param(object) for(int i=0;i<initial_params::num_initial_params;i++) {if(withinbound(0,(initial_params::varsptr[i])->phase_start, initial_params::current_phase)) { int sc= (initial_params::varsptr[i])->size_count(); if (sc>0) { write_input_log << "# " << initial_params::varsptr[i] ->label() << "\n" << object<<endl; } }}
	#define log_param(object) if (active(object)) write_input_log << "# " #object "\n" << object << endl;
  ofstream write_input_log("input.log");

 // void get_sel_changes(int& k);
  adstring_array stkname;
  adstring_array fshname;
  adstring_array indname;
  adstring truname;
  adstring simname;
  adstring model_name;
  adstring projfile_name;
  adstring datafile_name;
  adstring cntrlfile_name;
  adstring tmpstring;
  adstring repstring;
  adstring version_info;

FUNCTION dvariable get_spr_rates(double spr_percent,int istk)
  RETURN_ARRAYS_INCREMENT();
  dvar_matrix sel_tmp(1,nages,1,nfsh);
  sel_tmp.initialize();
  for (k=1;k<=nfsh;k++)
    if (sel_map(1,k) == istk)
      for (j=1;j<=nages;j++)
        sel_tmp(j,k) = sel_fsh(k,endyr,j); // NOTE uses last-year of fishery selectivity for projections.
  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
  {
    Fratio(k) = sum(F(k,endyr)) ;
    sumF(sel_map(1,k)) += Fratio(k) ;
  }
  for (k=1;k<=nfsh;k++)
    Fratio(k) /= sumF(sel_map(1,k));
  double df=1.e-3;
  dvariable F1 ;
  F1.initialize();
  F1 = .8*natmortprior(mort_map(istk,1)); //Ojo: first regime //nreg(istk)
  dvariable F2;
  dvariable F3;
  dvariable yld1;
  dvariable yld2;
  dvariable yld3;
  dvariable dyld;
  dvariable dyldp;
  // Newton Raphson stuff to go here
  for (int ii=1;ii<=6;ii++)
  {
    F2     = F1 + df;
    F3     = F1 - df;
    yld1   = -1000*square(log(spr_percent/spr_ratio(F1, sel_tmp,styr,istk)));
    yld2   = -1000*square(log(spr_percent/spr_ratio(F2, sel_tmp,styr,istk)));
    yld3   = -1000*square(log(spr_percent/spr_ratio(F3, sel_tmp,styr,istk)));
    dyld   = (yld2 - yld3)/(2*df);                          // First derivative (to find the root of this)
    dyldp  = (yld3-(2*yld1)+yld2)/(df*df);  // Newton-Raphson approximation for second derivitive
    F1    -= dyld/dyldp;
  }
  RETURN_ARRAYS_DECREMENT();
  return(F1);
  
FUNCTION dvariable spr_ratio(dvariable trial_F,dvar_matrix sel_tmp,int iyr,int istk)
  dvariable SBtmp;
  dvar_vector Ntmp(1,nages);
  dvar_vector srvtmp(1,nages);
  SBtmp.initialize();
  Ntmp.initialize();
  srvtmp.initialize();
  dvar_matrix Ftmp(1,nages,1,nfsh); // note that this is in reverse order of usual indexing (age, fshery)
  Ftmp = sel_tmp;
  for (j=1;j<=nages;j++) 
  {
    Ftmp(j) = elem_prod(Ftmp(j), trial_F * Fratio);
    srvtmp(j)  = mfexp(-sum(Ftmp(j)) - M(istk,iyr,j));
  }
  Ntmp(1)=1.;
  j=1;
  SBtmp  += Ntmp(j)*wt_mature(istk,j)*pow(srvtmp(j),spmo_frac);
  for (j=2;j<nages;j++)
  {
    Ntmp(j) = Ntmp(j-1)*srvtmp(j-1);
    SBtmp  += Ntmp(j)*wt_mature(istk,j)*pow(srvtmp(j),spmo_frac);
  }
  Ntmp(nages)=Ntmp(nages-1)*srvtmp(nages-1)/(1.-srvtmp(nages));
  SBtmp  += Ntmp(nages)*wt_mature(istk,nages)*pow(srvtmp(nages),spmo_frac);
  return(SBtmp/phizero(cum_regs(istk)+yy_sr(istk,iyr)));
  
FUNCTION dvariable spr_unfished(int istk,int i)
  dvariable Ntmp;
  dvariable SBtmp;
  SBtmp.initialize();
  Ntmp = 1.;
  for (j=1;j<nages;j++)
  {
    SBtmp += Ntmp*wt_mature(istk,j)*exp(-spmo_frac * M(istk,i,j));
    Ntmp  *= mfexp( -M(istk,i,j));
  }
  Ntmp    /= (1.-exp(-M(istk,i,nages)));
  SBtmp += Ntmp*wt_mature(istk,nages)*exp(-spmo_frac * M(istk,i,nages));
  return(SBtmp);

FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  dvar_vector sumF(1,nstk);
  sumF.initialize();
  for (k=1;k<=nfsh;k++)
  {
    Fratio(k) = sum(F(k,endyr)) ;
    sumF(sel_map(1,k)) += Fratio(k) ;
  }
  for (k=1;k<=nfsh;k++)
    Fratio(k) /= sumF(sel_map(1,k));

  for (s=1;s<=nstk;s++)
  {
    F35_est(s) = get_spr_rates(.35,s);
    F50_est(s) = get_spr_rates(.50,s);
    F40_est(s) = get_spr_rates(.40,s);
  }

  for (k=1;k<=nfsh;k++)
  {
    F50(k) = F50_est(sel_map(1,k)) * (Fratio(k));
    F40(k) = F40_est(sel_map(1,k)) * (Fratio(k));
    F35(k) = F35_est(sel_map(1,k)) * (Fratio(k));
  }
  cout << F50<<endl<<F40<<endl<<F35<<endl;

FUNCTION void writerep(dvariable& tmp,adstring& tmpstring)
  cout <<tmpstring<<endl<<endl;
  tmpstring = printf("3.5%f",value(tmp));

FUNCTION dvariable SolveF2(const int& iyr, const dvar_vector& N_tmp, const double&  TACin, const int& istk)
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc; 
  dvar_matrix Fratsel(1,nfsh,1,nages);
  dvar_vector M_tmp(1,nages) ;
  dvar_vector Z_tmp(1,nages) ;
  dvar_vector S_tmp(1,nages) ;
  dvar_vector Ftottmp(1,nages);
  dvariable btmp =  N_tmp * elem_prod(sel_fsh(1,iyr),wt_pop(istk));
  dvariable ftmp;
  M_tmp = M(istk,iyr);
  ftmp = TACin/btmp;
    for (k=1;k<=nfsh;k++)
      Fratsel(k) = Fratio(k)*sel_fsh(k,iyr);
    for (int ii=1;ii<=5;ii++)
    {
      Ftottmp.initialize();
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == istk)
          Ftottmp += ftmp*Fratsel(k);
  
      Z_tmp = Ftottmp  + M_tmp; 
      S_tmp = mfexp( -Z_tmp );
      cc = 0.0;
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == istk)
          cc += wt_fsh(k,endyr) * elem_prod(elem_div(ftmp*Fratsel(k),  Z_tmp),elem_prod(1.-S_tmp,N_tmp)); // Catch equation (vectors)
  
      dd = cc / TACin - 1.;
      if (dd<0.) dd *= -1.;
      ftmp += (TACin-cc) / btmp;
    }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp);

FUNCTION dvar_vector SolveF2(const int& iyr, const dvector&  Catch, const int& istk)
  // Returns vector of F's (given year) by fleet
  // Requires: N and fleet specific wts & selectivities at age, catch 
  // iterate to get Z's right
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc; 
  dvar_matrix  seltmp(1,nfsh,1,nages);
  dvar_matrix  wt_tmp(1,nfsh,1,nages);
  dvar_matrix Fratsel(1,nfsh,1,nages);
  dvar_vector N_tmp = natage(istk,iyr);
  dvar_vector M_tmp(1,nages) ;
  dvar_vector Z_tmp(1,nages) ;
  dvar_vector S_tmp(1,nages) ;
  dvar_vector Ftottmp(1,nages);
  dvar_vector Frat(1,nfsh);
  dvar_vector btmp(1,nfsh);
  dvar_vector ftmp(1,nfsh);
  dvar_vector hrate(1,nfsh);
  btmp.initialize(); 
  M_tmp = M(istk,iyr);
  // Initial guess for Fratio
  for (k=1;k<=nfsh;k++)
  {
    seltmp(k)= sel_fsh(k,iyr); // Selectivity
    wt_tmp(k)= wt_fsh(k,iyr); // 
    btmp(k)  =  N_tmp * elem_prod(seltmp(k),wt_tmp(k));
    hrate(k) = Catch(k)/btmp(k);
    Frat(k)  = Catch(k)/sum(Catch);
    Fratsel(k) = Frat(k)*seltmp(k);
    ftmp(k) = 1.1*(1.- posfun(1.-hrate(k),.10,fpen(4)));
  }
  // Initial fleet-specific F
  // iterate to balance effect of multiple fisheries...........
  for (int kk=1;kk<=nfsh;kk++) 
  {
    if (sel_map(1,kk) == istk)
    {
      for (k=1;k<=nfsh;k++)
      {
        if (sel_map(1,k) == istk)
        {
          if (hrate(k) <.9999) 
          {
            for (int ii=1;ii<=8;ii++)
            {
              Ftottmp.initialize();
              for (int aaa=1;aaa<=nages;aaa++)
                 for (int kkk=1;kkk<=nfsh;kkk++)
                   if (sel_map(1,kkk) == istk)
                     Ftottmp(aaa)   += ftmp(kkk)*Fratsel(kkk,aaa);
              Z_tmp     = Ftottmp  + M_tmp; 
              S_tmp     = mfexp( -Z_tmp );
              cc        = wt_tmp(k) * elem_prod(elem_div(ftmp(k)*Fratsel(k),  Z_tmp),elem_prod(1.-S_tmp,N_tmp)); // Catch equation (vectors)
              ftmp(k)  += ( Catch(k)-cc ) / btmp(k);
            }
            Frat(k)    = ftmp(k)/sum(ftmp);
            Fratsel(k) = Frat(k)*seltmp(k);
          }
        }
      }
    }
  }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp);

FUNCTION Write_SimDatafile
  {
  int nsims;
  // get the number of simulated datasets to create...
  ifstream sim_in("nsims.dat"); sim_in >> nsims; sim_in.close();
  ofstream SimDB("simout.dat",ios::app); 
  ofstream TruDB("truout.dat",ios::app); 
  // compute the autocorrelation term for residuals of fit to indices...
  // for (k=1;k<=nind;k++) ac(k) = get_AC(k);
  int nyrs_fsh_age_sim= 1+endyr-styr;
  int nyrs_ind_sim    = 1+endyr-styr;
  int nyrs_ind_age_sim= 1+endyr-styr;
  ivector yrs_fsh_age_sim(1,nyrs_fsh_age_sim);
  ivector yrs_ind_sim(1,nyrs_ind_sim);
  ivector yrs_ind_age_sim(1,nyrs_ind_age_sim);
  yrs_fsh_age_sim.fill_seqadd(1977,1);
  yrs_ind_sim.fill_seqadd(1977,1);
  yrs_ind_age_sim.fill_seqadd(1977,1);
  ivector n_sample_fsh_age_sim(1,nyrs_fsh_age_sim);
  ivector n_sample_ind_age_sim(1,nyrs_ind_age_sim);
  dvector new_ind_sim(1,nyrs_ind_sim);
  dmatrix sim_rec_devs(1,nstk,styr_rec,endyr);
  dmatrix sim_Sp_Biom(1,nstk,styr_rec,endyr);
  d3_array sim_natage(1,nstk,styr_rec,endyr,1,nages);
  d3_array catagetmp(1,nstk,styr,endyr,1,nages);
  dmatrix sim_catchbio(1,nstk,styr,endyr);
  dvector survtmp(1,nstk);
  for (s=1;s<=nstk;s++)
    survtmp(s) = value(mfexp(-natmort(s,styr)));
  Ftot.initialize();// Ojo
  for (k=1;k<=nfsh;k++) Ftot(sel_map(1,k)) += F(k);

  for (int isim=1;isim<=nsims;isim++)
  {
    new_ind_sim.initialize();
    sim_natage.initialize();
    // Start w/ simulated population
    // Simulate using new recruit series (same F's)
    // fill vector with unit normal RVs
    sim_rec_devs.fill_randn(rng);
    for (s=1;s<=nstk;s++)
      for (i=styr_rec;i<=endyr;i++)
        sim_rec_devs(s,i) *= value(sigmar(rec_map(s,yy_sr(s,i))));
    for (s=1;s<=nstk;s++)
    {
      sim_natage(s,styr_rec,1) = value(Rzero(cum_regs(s)+1))*exp(sim_rec_devs(s,styr_rec));
      for (j=2; j<=nages; j++)
        sim_natage(s,styr_rec,j) = sim_natage(s,styr_rec,j-1) * survtmp(s);
      sim_natage(s,styr_rec,nages) /= (1.-survtmp(s));
    }
  
    // Simulate population w/ process errors in recruits
    for (s=1;s<=nstk;s++)
    {
      for (i=styr_rec;i<=endyr;i++)
      {
        sim_Sp_Biom(s,i) = sim_natage(s,i)*pow(survtmp(s),spmo_frac) * wt_mature(s); 
        if (i>styr_rec+rec_age)
          sim_natage(s,i,1)          = value(SRecruit(sim_Sp_Biom(s,i-rec_age),cum_regs(s)+yy_sr(s,i)))*mfexp(sim_rec_devs(s,i)); 
        else
          sim_natage(s,i,1)          = value(SRecruit(sim_Sp_Biom(s,i),cum_regs(s)+yy_sr(s,i)))*mfexp(sim_rec_devs(s,i)); 
  
        if (i>=styr)
        {
          // apply estimated survival rates
          sim_Sp_Biom(s,i)          = value( elem_prod(sim_natage(s,i),pow(S(s,i),spmo_frac)) * wt_mature(s)); 
          catagetmp(s,i)            = value( elem_prod(elem_div(Ftot(s,i),Z(s,i)),elem_prod(1.-S(s,i),sim_natage(s,i))));
          sim_catchbio(s,i)         = catagetmp(s,i)*wt_fsh(1,i); //Ojo
          if (i<endyr)
          {
            sim_natage(s,i+1)(2,nages) = value( ++elem_prod(sim_natage(s,i)(1,nages-1),S(s,i)(1,nages-1)));
            sim_natage(s,i+1,nages)   += value( sim_natage(s,i,nages)*S(s,i,nages));
          }
        }
        else
        {
          if (i<endyr)
          {
            sim_natage(s,i+1)(2,nages) = ++(sim_natage(s,i)(1,nages-1) * survtmp(s));
            sim_natage(s,i+1,nages)   += sim_natage(s,i,nages)*survtmp(s);
          }
        }
      }
    }
  
    //===============================================
    //Now write from simulated population
    //
    // Create the name of the simulated dataset
    simname = "sim_"+ str(isim) + ".dat";
    truname = "tru_"+ str(isim) + ".dat";
    ofstream trudat(truname);
    truth(Rzero);
    truth(Fmsy);
    truth(MSY);
    dmatrix ntmp(1,nstk,1,nages);
    dmatrix seltmp(1,nfsh,1,nages);
    dmatrix Fatmp(1,nfsh,1,nages);
    dmatrix Ztmp(1,nstk,1,nages);
    seltmp.initialize();
    Fatmp.initialize();
    Ztmp.initialize();
    ntmp.initialize();
    for (k=1;k<=nfsh;k++)
      seltmp(k) = value(sel_fsh(k,endyr));
    for (s=1;s<=nstk;s++)
      Ztmp(s) = value(natmort(s,styr));
    for (k=1;k<=nfsh;k++)
    { 
      Fatmp(k) = value(Fratio(k) * Fmsy(sel_map(1,k)) * seltmp(k));
      Ztmp(sel_map(1,k))    += Fatmp(k);
    } 
    dmatrix survmsy = exp(-Ztmp);
    for (s=1;s<=nstk;s++)
    {
      ntmp(s,1) = value(Rmsy(s));
      for (j=2;j<=nages;j++) 
        ntmp(s,j) = ntmp(s,j-1)*survmsy(s,j-1);
      ntmp(s,nages) /= (1-survmsy(s,nages));
    }
    // dvariable phi    = elem_prod( ntmp , pow(survmsy,spmo_frac ) )* wt_mature;
    truth(Rmsy);
    truth(seltmp);
    dvector SurvBmsy(1,nstk);
    double q_ind_sim=value(mean(q_ind(1))); //Ojo
    for (s=1;s<=nstk;s++)
      SurvBmsy(s) = value(elem_prod(wt_ind(1,endyr),elem_prod(pow(survmsy(s),ind_month_frac(1)), ntmp(s))) * q_ind_sim*sel_ind(1,endyr)); //Ojo
    truth(ntmp);
    dvector Cmsy(1,nstk);
    for (s=1;s<=nstk;s++)
      Cmsy(s) = value(yield(Fratio, Fmsy(s), s));
    truth(Cmsy);
    // Now do OFL for next year...
    for (s=1;s<=nstk;s++)
    {
      ntmp(s,1)        = value(SRecruit(sim_Sp_Biom(s,endyr+1-rec_age),cum_regs(s)+yy_sr(s,endyr+1)));
      ntmp(s)(2,nages) = value( ++elem_prod(sim_natage(s,endyr)(1,nages-1),S(s,endyr)(1,nages-1)));
      ntmp(s,nages)   += value( sim_natage(s,endyr,nages)*S(s,endyr,nages));
    }
    dvector ctmp(1,nages);
    ctmp.initialize();
    OFL.initialize();
    for (k=1;k<=nfsh;k++)
    {
      for ( j=1 ; j <= nages; j++ )
        ctmp(j)      = ntmp(sel_map(1,k),j) * Fatmp(k,j) * (1. - survmsy(sel_map(1,k),j)) / Ztmp(sel_map(1,k),j);
      OFL(sel_map(1,k))  += wt_fsh(k,endyr) * ctmp;
    }
    dvector NextSurv(1,nstk);
    dvector NextSSB(1,nstk);
    for (s=1;s<=nstk;s++)
    {
      NextSurv(s) = value(elem_prod(wt_ind(1,endyr),elem_prod(pow(survmsy(s),ind_month_frac(1)), ntmp(s))) * 
                   q_ind_sim*sel_ind(1,endyr)); //Ojo
      NextSSB(s)  = elem_prod(ntmp(s), pow(survmsy(s),spmo_frac)) * wt_mature(s); 
    }
    // Catch at following year for Fmsy
    truth(OFL);
    truth(SurvBmsy);
    truth(steepness);
    truth(natmort);
    truth(sim_natage);
    truth(sim_Sp_Biom);
    // Open the simulated dataset for writing
    ofstream simdat(simname);
    simdat << "# first year" <<endl;
    simdat << styr <<endl;
    simdat << "# Last  year" <<endl;
    simdat << endyr <<endl;
    simdat << "# age recruit" <<endl;
    simdat << rec_age <<endl;
    simdat << "# oldest age" <<endl;
    simdat << oldest_age <<endl;
    simdat << "# Number of stocks " <<endl;
    simdat << nstk <<endl;                                   
    simdat << stknameread <<endl;
    simdat << "# Number of fisheries " <<endl;
    simdat << nfsh <<endl;                                   
    simdat << fshnameread <<endl;                                   
    simdat << "# Catch biomass by stock " <<endl;
    for (s=1;s<=nstk;s++)
    {
      simdat << "# " <<stkname(s) <<" " << s <<endl;
      simdat << sim_catchbio(s) <<endl;
    }
    simdat << "# Catch biomass uncertainty by fishery (std errors)" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k) <<" " << k <<endl;
      simdat << catch_bio_sd(k) <<endl;   
    }
    simdat << "# number of years for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << nyrs_fsh_age_sim <<endl;
    }
    simdat << "# years for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << yrs_fsh_age_sim  <<endl;
    }
    simdat << "# sample sizes for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      n_sample_fsh_age_sim = mean(n_sample_fsh_age(k));
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << n_sample_fsh_age_sim         <<endl;    
    }
    simdat << "# Observed age compositions for fishery" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      dvector p(1,nages);
      double Ctmp; // total catch
      dvector freq(1,nages);
      simdat << "# " << fshname(k) <<endl;
      for (i=1;i<=nyrs_fsh_age_sim;i++)
      {
        int iyr = yrs_fsh_age_sim(i);
        // Add noise here
        freq.initialize();
        ivector bin(1,n_sample_fsh_age_sim(i));
        p  = catagetmp(sel_map(1,k),iyr);
        p /= sum(p);
        bin.fill_multinomial(rng,p); // fill a vector v
        for (int j=1;j<=n_sample_fsh_age_sim(i);j++)
          freq(bin(j))++;
        // Apply ageing error to samples..............
        // p = age_err *freq/sum(freq); 
        p = freq/sum(freq); 
        simdat << p  <<endl;
        // Compute total catch given this sample size for catch-age
        Ctmp = sim_catchbio(sel_map(1,k),iyr) / (p*wt_fsh(k,iyr)); 
        // Simulated catage = proportion sampled
        // sim_catage(k,i) = p * Ctmp;
      }
    }
    simdat << "# Annual wt-at-age for fishery" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << (k) <<endl;
      // Add noise here
      simdat << wt_fsh(k)  <<endl;  
    }
    simdat << "# number of indices" <<endl;
    simdat << nind <<endl;                                   
    simdat << indnameread <<endl;                                   
    simdat << "# Number of years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << nyrs_ind_sim  <<endl;                   
    }
    simdat << "# Years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << yrs_ind_sim <<endl;         
    }
    simdat << "# Month that index occurs "<<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << mo_ind(k) <<endl;
    }
    simdat << "# values for indices (annual)"<<endl;
    // note assumes only one index...
    double ind_sigma;
    dvector ind_devs(1,nyrs_ind_sim);
    for (k=1;k<=nind;k++)
    {
      ind_sigma = mean(obs_lse_ind(k)) ;
      ind_sigma = 0.10 ;
      simdat << "# " <<indname(k)<< " " << k <<endl;
      // Add noise here
      // fill vector with unit normal RVs
      ind_devs.fill_randn(rng);
      ind_devs *= ind_sigma ;
      for (i=1;i<=nyrs_ind_sim;i++)
      {
        int iyr=yrs_ind_sim(i);
        //uncorrelated...corr_dev(k,i) = ac(k) * corr_dev(k,i-1) + sqrt(1.-square(ac(k))) * corr_dev(k,i);
        new_ind_sim(i) = mfexp(ind_devs(i) - ind_sigma/2.) * value(elem_prod(wt_ind(k,iyr),elem_prod(pow(S(sel_map(1,k+nfsh),iyr),ind_month_frac(k)), 
                        sim_natage(sel_map(1,k+nfsh),iyr))) * q_ind_sim*sel_ind(k,iyr)); 
      }
      simdat << new_ind_sim     <<endl;
      dvector ExactSurvey = elem_div(new_ind_sim,exp(ind_devs-ind_sigma/2.));
      truth(ExactSurvey);
    }
    simdat << "# standard errors for indices (by year) " <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< " " << k <<endl;
      // simdat << new_ind_sim*mean(elem_div(obs_se_ind(k),obs_ind(k)))  <<endl;
      simdat << new_ind_sim*ind_sigma  <<endl; //Ojo
    }
    simdat << "# Number of years of age data available for index" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< " " << k <<endl;
      simdat << nyrs_ind_age_sim <<endl;
    }
    simdat << "# Years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      simdat << yrs_ind_age_sim <<endl;
    }
    simdat << "# Sample sizes for age data from indices" <<endl;
    for (k=1;k<=nind;k++)
    {
      n_sample_ind_age_sim = mean(n_sample_ind_age(k));
      simdat << "# " <<indname(k)<< endl;
      simdat << n_sample_ind_age_sim <<endl;
    }
    simdat << "# values of proportions at age in index" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      dvector p(1,nages);
      dvector freq(1,nages);
      for (i=1;i<=nyrs_ind_age_sim;i++)
      {
        int iyr = yrs_ind_age_sim(i);
        // Add noise here
        freq.initialize();
        ivector bin(1,n_sample_ind_age_sim(i));
        // p = age_err * value(elem_prod( elem_prod(pow(S(iyr),ind_month_frac(k)), sim_natage(iyr))*q_ind_sim , sel_ind(k,iyr))); 
        p = value(elem_prod( elem_prod(pow(S(sel_map(1,k+nfsh),iyr),ind_month_frac(k)), sim_natage(sel_map(1,k+nfsh),iyr))*q_ind_sim , sel_ind(k,iyr))); 
        p /= sum(p);
        // fill vector with multinomial samples
        bin.fill_multinomial(rng,p); // fill a vector v
        for (int j=1;j<=n_sample_ind_age_sim(i);j++)
          freq(bin(j))++;
        simdat << "# " <<indname(k)<< " year: "<< iyr<< endl;
        simdat << freq/sum(freq) <<endl;
      }
    }
    simdat << "# Mean wts at age for indices" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      // Could add noise here
      simdat <<  wt_ind(k)  <<endl;
    }
  
    simdat << "# Population mean wt at age" <<endl;
    simdat << wt_pop <<endl;
  
    simdat << "# Population maturity at age" <<endl;
    simdat << maturity  <<endl;
  
    simdat << "# Peak spawning month" <<endl;
    simdat << spawnmo <<endl;
  
    simdat << "# ageing error " <<endl;
    simdat << age_err <<endl;

    simdat <<endl<<endl<<"Additional output"<<endl;
    simdat << "# Fishery_Effort " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      dvector ran_fsh_vect(styr,endyr);
      // fill vector with unit normal RVs
      ran_fsh_vect.fill_randn(rng);
      // Sigma on effort is ~15% white noise (add red noise later)
      ran_fsh_vect *= 0.15; 
      dvector avail_biom(styr,endyr);
      for (i=styr;i<=endyr;i++)
      {
        avail_biom(i) = wt_fsh(k,i)*value(elem_prod(sim_natage(sel_map(1,k),i),sel_fsh(k,i))); 
      }
      act_eff(k) = elem_prod(exp(ran_fsh_vect), (elem_div(catch_bio(k), avail_biom)) );
      // Normalize effort
      act_eff(k) /= mean(act_eff(k));
      for (i=styr;i<=endyr;i++)
        simdat<<fshname(k)<<" "<<i<<" "<<act_eff(k,i) <<endl;
    }
    simdat << "# Fishery catch-at-age " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << "Fishery Year "<<age_vector << endl;
      for (i=1;i<=nyrs_fsh_age(k);i++)
        simdat<<fshname(k)<<" "<<yrs_fsh_age(k,i)<<" "<<catagetmp(sel_map(1,k),yrs_fsh_age(k,i)) <<endl;
    }
    // Write simple file by simulation
    dvector ExactSurvey = elem_div(new_ind_sim,exp(ind_devs-ind_sigma/2.));
    for (s=1;s<=nstk;s++)
    {
      for (i=styr;i<=endyr;i++)
      {
        SimDB<<model_name<<" "<<isim<<" "<< i<<" "<<
          sim_catchbio(s,i)     <<" "<< 
          new_ind_sim(i-styr+1) <<" "<< 
          new_ind_sim(i-styr+1)*ind_sigma  <<endl;
        TruDB<<model_name<<" " <<isim<<" "<< i<<" "<<
          sim_catchbio(s,i)    <<" "<< 
          sim_natage(s,i,1)    <<" "<< 
          sim_Sp_Biom(s,i)     <<" "<< 
          ExactSurvey(i-styr+1)<<" "<<
          steepness(cum_regs(s)+yy_sr(s,i))<<" "<< 
          Bmsy(s)              <<" "<< 
          MSYL(s)              <<" "<< 
          MSY(s)               <<" "<< 
          SurvBmsy(s)          <<" "<<
          endl;
      }
      TruDB<<model_name<<" "<<isim<<" "<< endyr+1<<" "<<
          OFL(s)               <<" "<< 
          SRecruit(sim_Sp_Biom(s,endyr+1-rec_age),cum_regs(s)+yy_sr(s,endyr+1))<<" "<<
          sim_Sp_Biom(s,endyr) <<" "<< 
          NextSurv             <<" "<< 
          steepness(cum_regs(s)+yy_sr(s,i))<<" "<< 
          Bmsy(s)              <<" "<< 
          MSYL(s)              <<" "<< 
          MSY(s)               <<" "<< 
          SurvBmsy(s)          <<" "<<
          endl;
    }

    trudat.close();
  }
  SimDB.close();
  TruDB.close();
  exit(1);
  // End of simulating datasets...................
  }

FUNCTION Write_Datafile
  dmatrix new_ind(1,nind,1,nyrs_ind);
  new_ind.initialize();
  int nsims;
  // get the number of simulated datasets to create...
  ifstream sim_in("nsims.dat"); sim_in >> nsims; sim_in.close();
  // compute the autocorrelation term for residuals of fit to indices...
  for (k=1;k<=nind;k++)
    ac(k) = get_AC(k);
  for (int isim=1;isim<=nsims;isim++)
  {
    // Create the name of the simulated dataset
    simname = "sim_"+ str(isim) + ".dat";
    // Open the simulated dataset for writing
    ofstream simdat(simname);
    simdat << "# first year" <<endl;
    simdat << styr <<endl;
    simdat << "# Last  year" <<endl;
    simdat << endyr <<endl;
    simdat << "# age recruit" <<endl;
    simdat << rec_age <<endl;
    simdat << "# oldest age" <<endl;
    simdat << oldest_age <<endl;
    simdat << "# Number of fisheries " <<endl;
    simdat << nfsh <<endl;                                   
    simdat << fshnameread <<endl;                                   
    simdat << "# Catch biomass by fishery " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k) <<" " << k <<endl;
      simdat << catch_bio(k) <<endl;
    }
    simdat << "# Catch biomass uncertainty by fishery (std errors)" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k) <<" " << k <<endl;
      simdat << catch_bio_sd(k) <<endl;   
    }
    simdat << "# number of years for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << nyrs_fsh_age(k) <<endl;
    }
    simdat << "# years for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << yrs_fsh_age(k)  <<endl;
    }
    simdat << "# sample sizes for fishery age data " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << n_sample_fsh_age(k)  <<endl;    
    }
    simdat << "# Observed age compositions for fishery" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      dvector p(1,nages);
      double Ctmp; // total catch
      dvector freq(1,nages);
      simdat << "# " << fshname(k) <<endl;
      for (i=1;i<=nyrs_fsh_age(k);i++)
      {
        int iyr = yrs_fsh_age(k,i);
        // Add noise here
        freq.initialize();
        ivector bin(1,n_sample_fsh_age(k,i));
        p  = value(catage(k,iyr));
        p /= sum(p);
        bin.fill_multinomial(rng,p); // fill a vector v
        for (int j=1;j<=n_sample_fsh_age(k,i);j++)
          freq(bin(j))++;
        // Apply ageing error to samples..............
        p = age_err *freq/sum(freq); 
        simdat << p  <<endl;
        // Compute total catch given this sample size
        Ctmp = catch_bio(k,iyr) / (p*wt_fsh(k,iyr)); 
        // Simulated catage = proportion sampled
        catage(k,i) = p * Ctmp;
      }
    }
    simdat << "# Annual wt-at-age for fishery" <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << (k) <<endl;
      // Add noise here
      simdat << wt_fsh(k)  <<endl;  
    }
    simdat << "# number of indices" <<endl;
    simdat << nind <<endl;                                   
    simdat << indnameread <<endl;                                   
    simdat << "# Number of years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << nyrs_ind(k)  <<endl;                   
    }
    simdat << "# Years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << yrs_ind(k)  <<endl;         
    }
    simdat << "# Month that index occurs "<<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " << indname(k) <<endl;
      simdat << mo_ind(k) <<endl;
    }
    simdat << "# values for indices (annual)"<<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< " " << k <<endl;
      // Add noise here
      dvector ran_ind_vect(1,nyrs_ind(k));
      // fill vector with unit normal RVs
      ran_ind_vect.fill_randn(rng);
      // do first year uncorrelated
      i=1;
      int iyr=yrs_ind(k,i);
      corr_dev(k)  = ran_ind_vect;
      new_ind(k,i) = mfexp(corr_dev(k,i) * obs_lse_ind(k,i) ) * 
                     value(elem_prod(wt_ind(k,iyr),elem_prod(pow(S(sel_map(1,k+nfsh),iyr),ind_month_frac(k)), natage(sel_map(1,k+nfsh),iyr)))*
                     q_ind(k,i)*sel_ind(k,iyr)); 
      // do next years correlated with previous
      for (i=2;i<=nyrs_ind(k);i++)
      {
        iyr=yrs_ind(k,i);
        corr_dev(k,i) = ac(k) * corr_dev(k,i-1) + sqrt(1.-square(ac(k))) * corr_dev(k,i);
        new_ind(k,i) = mfexp(corr_dev(k,i) * obs_lse_ind(k,i) ) * 
                        value(elem_prod(wt_ind(k,iyr),elem_prod(pow(S(sel_map(1,k+nfsh),iyr),ind_month_frac(k)), 
                        natage(sel_map(1,k+nfsh),iyr))) * q_ind(k,i)*sel_ind(k,iyr)); 
      }
      simdat << new_ind(k)      <<endl;
    }
    simdat << "# standard errors for indices (by year) " <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< " " << k <<endl;
      simdat << obs_se_ind(k)  <<endl;
    }
    simdat << "# Number of years of age data available for index" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< " " << k <<endl;
      simdat << nyrs_ind_age(k)  <<endl;
    }
    simdat << "# Years of index values (annual)" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      simdat << yrs_ind_age(k)  <<endl;
    }
    simdat << "# Sample sizes for age data from indices" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      simdat << n_sample_ind_age(k)  <<endl;
    }
    simdat << "# values of proportions at age in index" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      dvector p(1,nages);
      dvector freq(1,nages);
      for (i=1;i<=nyrs_ind_age(k);i++)
      {
        int iyr = yrs_ind_age(k,i);
        // Add noise here
        freq.initialize();
        ivector bin(1,n_sample_ind_age(k,i));
        p = age_err * value(elem_prod( elem_prod(pow(S(sel_map(1,k+nfsh),iyr),ind_month_frac(k)), natage(sel_map(1,k+nfsh),iyr))*q_ind(k,i) , sel_ind(k,iyr))); 
        p /= sum(p);
        // fill vector with multinomial samples
        bin.fill_multinomial(rng,p); // fill a vector v
        for (int j=1;j<=n_sample_ind_age(k,i);j++)
          freq(bin(j))++;
        simdat << "# " <<indname(k)<< " year: "<< iyr<< endl;
        simdat << freq/sum(freq) <<endl;
      }
    }
    simdat << "# Mean wts at age for indices" <<endl;
    for (k=1;k<=nind;k++)
    {
      simdat << "# " <<indname(k)<< endl;
      // Could add noise here
      simdat <<  wt_ind(k)  <<endl;
    }
  
    simdat << "# Population mean wt at age" <<endl;
    simdat << wt_pop <<endl;
  
    simdat << "# Population maturity at age" <<endl;
    simdat << maturity  <<endl;
  
    simdat << "# Peak spawning month" <<endl;
    simdat << spawnmo <<endl;
  
    simdat << "# ageing error " <<endl;
    simdat << age_err <<endl;

    simdat <<endl<<endl<<"Additional output"<<endl;
    simdat << "# Fishery_Effort " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      dvector ran_fsh_vect(styr,endyr);
      // fill vector with unit normal RVs
      ran_fsh_vect.fill_randn(rng);
      // Sigma on effort is ~15% white noise (add red noise later)
      ran_fsh_vect *= 0.15; 
      dvector avail_biom(styr,endyr);
      for (i=styr;i<=endyr;i++)
      {
        avail_biom(i) = wt_fsh(k,i)*value(elem_prod(natage(sel_map(1,k),i),sel_fsh(k,i))); 
      }
      act_eff(k) = elem_prod(exp(ran_fsh_vect), (elem_div(catch_bio(k), avail_biom)) );
      // Normalize effort
      act_eff(k) /= mean(act_eff(k));
      for (i=styr;i<=endyr;i++)
        simdat<<fshname(k)<<" "<<i<<" "<<act_eff(k,i) <<endl;
    }
    simdat << "# Fishery catch-at-age " <<endl;
    for (k=1;k<=nfsh;k++)
    {
      simdat << "# " <<fshname(k)<< " " << k <<endl;
      simdat << "Fishery Year "<<age_vector << endl;
      for (i=1;i<=nyrs_fsh_age(k);i++)
        simdat<<fshname(k)<<" "<<yrs_fsh_age(k,i)<<" "<<catage(k,i) <<endl;
    }
  }
  exit(1);
  // End of simulating datasets...................
  

FUNCTION Write_R
  adstring report_name;
  for (s=1;s<=nstk;s++)
  {
    report_name = "For_R_"+ str(s) + ".rep";
    ofstream R_report(report_name);
    R_report<< "$repl_yld"<<endl<<repl_yld(s)<<endl; 
    R_report<< "$repl_SSB"<<endl<<repl_SSB(s)<<endl; 
    for (int k=1;k<=nfsh;k++)
    {
      if (nyrs_fsh_age(k)>0)
        R_report<< "$FW_"<<fshname(k) << endl<< Francis_wts(k)<<endl; 
    }
    for (int k=1;k<=nind;k++)
    {
      if (nyrs_ind_age(k)>0)
        R_report<< "$FW_"<<indname(k) << endl<< Francis_wts(k+nfsh)<<endl; 
    }
    R_report<<"$M"<<endl; 
    R_report<<M(s)<<endl;
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report<<"$q_"<<k<<endl; 
        for (i=1;i<nyrs_ind(k);i++)
        {        
          int iyr=yrs_ind(k,i);
          for (int ii=iyr;ii<yrs_ind(k,i+1);ii++)
            R_report<<ii<<" "<<pow(q_ind(k,i),q_power_ind(k))<<endl;
        }
        R_report<<yrs_ind(k,nyrs_ind(k))<<" "<<pow(q_ind(k,nyrs_ind(k)),q_power_ind(k))<<endl;
      }
    }
    R_report<<"$SurvNextYr"<<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report<< pred_ind_nextyr(k) <<" ";
      R_report<<endl;
    R_report<<"$Yr"<<endl; for (i=styr;i<=endyr;i++) R_report<<i<<" "; R_report<<endl;
    for (r=1;r<=nreg(s);r++)
    {
      R_report<<"$P_age2len_"<<growth_map(s,r)<<endl;
      for (j=1;j<=nages;j++)
        R_report<<P_age2len(growth_map(s,r),j)<<endl;
    }
    R_Report(len_bins);
    R_report<<"$TotF"<<endl << Ftot(s)<<endl;
    R_report<<"$TotBiom_NoFish"<<endl; for (i=styr;i<=endyr;i++) 
    {
      double lb=value(totbiom_NoFish(s,i)/exp(2.*sqrt(log(1+square(totbiom_NoFish.sd(s,i))/square(totbiom_NoFish(s,i))))));
      double ub=value(totbiom_NoFish(s,i)*exp(2.*sqrt(log(1+square(totbiom_NoFish.sd(s,i))/square(totbiom_NoFish(s,i))))));
      R_report<<i<<" "<<totbiom_NoFish(s,i)<<" "<<totbiom_NoFish.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report<<"$SSB_NoFishR"<<endl; for (i=styr+1;i<=endyr;i++) 
    {
      double lb=value(Sp_Biom_NoFishRatio(s,i)/exp(2.*sqrt(log(1+square(Sp_Biom_NoFishRatio.sd(s,i))/square(Sp_Biom_NoFishRatio(s,i))))));
      double ub=value(Sp_Biom_NoFishRatio(s,i)*exp(2.*sqrt(log(1+square(Sp_Biom_NoFishRatio.sd(s,i))/square(Sp_Biom_NoFishRatio(s,i))))));
      R_report<<i<<" "<<Sp_Biom_NoFishRatio(s,i)<<" "<< Sp_Biom_NoFishRatio.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }

    R_report<<"$TotBiom"<<endl; 
    for (i=styr;i<=endyr;i++) 
    {
      double lb=value(totbiom(s,i)/exp(2.*sqrt(log(1+square(totbiom.sd(s,i))/square(totbiom(s,i))))));
      double ub=value(totbiom(s,i)*exp(2.*sqrt(log(1+square(totbiom.sd(s,i))/square(totbiom(s,i))))));
      R_report<<i<<" "<<totbiom(s,i)<<" "<<totbiom.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }

    /*for (k=1;k<=5;k++){
      R_report<<"$SSB_fut_"<<k<<endl; 
      for (i=styr_fut;i<=endyr_fut;i++) 
      {
        double lb=value(SSB_fut(s,k,i)/exp(2.*sqrt(log(1+square(SSB_fut.sd(s,k,i))/square(SSB_fut(s,k,i))))));
        double ub=value(SSB_fut(s,k,i)*exp(2.*sqrt(log(1+square(SSB_fut.sd(s,k,i))/square(SSB_fut(s,k,i))))));
        R_report<<i<<" "<<SSB_fut(s,k,i)<<" "<<SSB_fut.sd(s,k,i)<<" "<<lb<<" "<<ub<<endl;
      }
    }*/
    R_report<<"$SSB_fut_1"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_1(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_1.sd(s,i))/square(SSB_fut_1(s,i))))));
      double ub=value(SSB_fut_1(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_1.sd(s,i))/square(SSB_fut_1(s,i))))));
      R_report<<i<<" "<<SSB_fut_1(s,i)<<" "<<SSB_fut_1.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report<<"$SSB_fut_2"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_2(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_2.sd(s,i))/square(SSB_fut_2(s,i))))));
      double ub=value(SSB_fut_2(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_2.sd(s,i))/square(SSB_fut_2(s,i))))));
      R_report<<i<<" "<<SSB_fut_2(s,i)<<" "<<SSB_fut_2.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report<<"$SSB_fut_3"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_3(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_3.sd(s,i))/square(SSB_fut_3(s,i))))));
      double ub=value(SSB_fut_3(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_3.sd(s,i))/square(SSB_fut_3(s,i))))));
      R_report<<i<<" "<<SSB_fut_3(s,i)<<" "<<SSB_fut_3.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report<<"$SSB_fut_4"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_4(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_4.sd(s,i))/square(SSB_fut_4(s,i))))));
      double ub=value(SSB_fut_4(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_4.sd(s,i))/square(SSB_fut_4(s,i))))));
      R_report<<i<<" "<<SSB_fut_4(s,i)<<" "<<SSB_fut_4.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report<<"$SSB_fut_5"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_5(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_5.sd(s,i))/square(SSB_fut_5(s,i))))));
      double ub=value(SSB_fut_5(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_5.sd(s,i))/square(SSB_fut_5(s,i))))));
      R_report<<i<<" "<<SSB_fut_5(s,i)<<" "<<SSB_fut_5.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
		/*
    R_report<<"$SSB_fut_6"<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut_6(s,i)/exp(2.*sqrt(log(1+square(SSB_fut_6.sd(s,i))/square(SSB_fut_6(s,i))))));
      double ub=value(SSB_fut_6(s,i)*exp(2.*sqrt(log(1+square(SSB_fut_6.sd(s,i))/square(SSB_fut_6(s,i))))));
      R_report<<i<<" "<<SSB_fut_6(s,i)<<" "<<SSB_fut_6.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
		*/

    double ctmp;
    for (k=1;k<=5;k++){
      R_report<<"$Catch_fut_"<<k<<endl; 
      for (i=styr_fut;i<=endyr_fut;i++) 
      {
        if (k==5) ctmp=0.;else ctmp=value(catch_future(s,k,i));
        R_report<<i<<" "<<ctmp<<endl;
      }
    }

    R_report<<"$SSB"<<endl; for (i=styr_sp;i<=endyr+1;i++) 
    {
      double lb=value(Sp_Biom(s,i)/exp(2.*sqrt(log(1+square(Sp_Biom.sd(s,i))/square(Sp_Biom(s,i))))));
      double ub=value(Sp_Biom(s,i)*exp(2.*sqrt(log(1+square(Sp_Biom.sd(s,i))/square(Sp_Biom(s,i))))));
      R_report<<i<<" "<<Sp_Biom(s,i)<<" "<<Sp_Biom.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }

    R_report<<"$R"<<endl; for (i=styr;i<=endyr;i++) 
    {
      double lb=value(recruits(s,i)/exp(2.*sqrt(log(1+square(recruits.sd(s,i))/square(recruits(s,i))))));
      double ub=value(recruits(s,i)*exp(2.*sqrt(log(1+square(recruits.sd(s,i))/square(recruits(s,i))))));
      R_report<<i<<" "<<recruits(s,i)<<" "<<recruits.sd(s,i)<<" "<<lb<<" "<<ub<<endl;
    }
    R_report << "$N"<<endl;
    for (i=styr;i<=endyr;i++) 
      R_report <<   i << " "<< natage(s,i) << endl;
      R_report   << endl;

    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report << "$F_age_"<< (k) <<""<< endl ;
        for (i=styr;i<=endyr;i++) 
          R_report <<i<<" "<<F(k,i)<<" "<< endl;
          R_report   << endl;
      }
    }

    R_report <<endl<< "$Fshry_names"<< endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report << fshname(k) << endl ;

    R_report <<endl<< "$Index_names"<< endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report << indname(k) << endl ;

    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        int ii=1;
        R_report <<endl<< "$Obs_Survey_"<< k <<""<< endl ;
        for (i=styr;i<=endyr;i++)
        {
          if (ii<=yrs_ind(k).indexmax())
          {
            if (yrs_ind(k,ii)==i)
            {
              double PearsResid   =  value((obs_ind(k,ii)-pred_ind(k,ii))/obs_se_ind(k,ii) );
              double lnPearsResid =  value((log(obs_ind(k,ii))-log(pred_ind(k,ii)))/obs_lse_ind(k,ii) );
              R_report << i<< " "<< obs_ind(k,ii)   <<" "<< 
                                    pred_ind(k,ii)  <<" "<< 
                                    obs_se_ind(k,ii)<<" "<<  
                                    PearsResid      <<" "<<
                                    lnPearsResid    << endl; //values of survey index value (annual)
              ii++;
            }
            // else
              // R_report << i<< " -1 "<< " "<< pred_ind(k,i)<<" -1 "<<endl;
          }
          // else
            // R_report << i<< " -1 "<< " "<< pred_ind(k,i)<<" -1 "<<endl;
        }
        R_report   << endl;
        R_report << endl<< "$Index_Q_"<<k<<endl;
        R_report<< q_ind(k) << endl;
      }
    }
    // R_report <<" SDNR1 "<< wt_srv1*std_dev(elem_div((pred_srv1(yrs_srv1)-obs_srv1_biom),obs_srv1_se))<<endl;
    R_report   << endl;
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        if (nyrs_fsh_age(k)>0) 
        { 
          R_report << "$pobs_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_age(k);i++) 
            R_report << yrs_fsh_age(k,i)<< " "<< oac_fsh(k,i) << endl;
          R_report   << endl;

          R_report << "$phat_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_age(k);i++) 
            R_report << yrs_fsh_age(k,i)<< " "<< eac_fsh(k,i) << endl;
            R_report   << endl;

          R_report << "$sdnr_age_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_age(k);i++) 
            R_report << yrs_fsh_age(k,i)<< " "<< sdnr( eac_fsh(k,i),oac_fsh(k,i),n_sample_fsh_age(k,i)) << endl;
          R_report   << endl;
        }
        if (nyrs_fsh_length(k)>0) 
        { 
          R_report << "$pobs_len_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_length(k);i++) 
            R_report << yrs_fsh_length(k,i)<< " "<< olc_fsh(k,i) << endl;
          R_report   << endl;

          R_report << "$phat_len_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_length(k);i++) 
            R_report << yrs_fsh_length(k,i)<< " "<< elc_fsh(k,i) << endl;
          R_report   << endl;

          R_report << "$sdnr_length_fsh_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_fsh_length(k);i++) 
            R_report << yrs_fsh_length(k,i)<< " "<< sdnr( elc_fsh(k,i),olc_fsh(k,i),n_sample_fsh_length(k,i)) << endl;
          R_report   << endl;
        }
      }
    }
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        if (nyrs_ind_age(k)>0) 
        { 
          R_report << "$pobs_ind_"<<(k)<<""<<  endl;
          for (i=1;i<=nyrs_ind_age(k);i++) 
            R_report << yrs_ind_age(k,i)<< " "<< oac_ind(k,i) << endl;
          R_report   << endl;
        
          R_report << "$phat_ind_"<<(k)<<""<<  endl;
          for (i=1;i<=nyrs_ind_age(k);i++) 
            R_report << yrs_ind_age(k,i)<< " "<< eac_ind(k,i) << endl;
          R_report   << endl;

          R_report << "$sdnr_age_ind_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_ind_age(k);i++) 
            R_report << yrs_ind_age(k,i)<< " "<< sdnr( eac_ind(k,i),oac_ind(k,i),n_sample_ind_age(k,i)) << endl;
          R_report   << endl;
        }
        if (nyrs_ind_length(k)>0) 
        { 
          R_report << "$pobs_len_ind_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_ind_length(k);i++) 
            R_report << yrs_ind_length(k,i)<< " "<< olc_ind(k,i) << endl;
          R_report   << endl;
          R_report << "$phat_len_ind_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_ind_length(k);i++) 
            R_report << yrs_ind_length(k,i)<< " "<< elc_ind(k,i) << endl;
          R_report   << endl;
          R_report << "$sdnr_length_ind_"<< (k) <<""<< endl;
          for (i=1;i<=nyrs_ind_length(k);i++) 
            R_report << yrs_ind_length(k,i)<< " "<< sdnr( elc_ind(k,i),olc_ind(k,i),n_sample_ind_length(k,i)) << endl;
          R_report   << endl;

        }
      }
    }
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report << endl<< "$Obs_catch_"<<(k) << endl;
        R_report << catch_bio(k) << endl;
        R_report   << endl;
        R_report << "$Pred_catch_" <<(k) << endl;
        R_report << pred_catch(k) << endl;
        R_report   << endl;
      }
    }

    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report << "$F_fsh_"<<(k)<<" "<<endl;
        for (i=styr;i<=endyr;i++)
        {
          R_report<< i<< " ";
          R_report<< mean(F(k,i)) <<" "<< mean(F(k,i))*max(sel_fsh(k,i)) << " ";
          R_report<< endl;
        }
      }
    }

    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report << endl<< "$sel_fsh_"<<(k)<<"" << endl;
        for (i=styr;i<=endyr;i++)
          R_report << k <<"  "<< i<<" "<<sel_fsh(k,i) << endl; 
        R_report   << endl;
      }
    }

    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report << endl<< "$sel_ind_"<<(k)<<"" << endl;
        for (i=styr;i<=endyr;i++)
          R_report << k <<"  "<< i<<" "<<sel_ind(k,i) << endl;
          R_report << endl;
      }
    }
    R_report << endl<< "$Stock_Rec"<< endl;
    for (i=styr_rec;i<=endyr;i++)
      if (active(log_Rzero(cum_regs(s)+yy_sr(s,i))))
        R_report << i<< " "<<Sp_Biom(s,i-rec_age)<< " "<< SRecruit(Sp_Biom(s,i-rec_age),cum_regs(s)+yy_sr(s,i))<< " "<< mod_rec(s,i)<<endl;
      else 
        R_report << i<< " "<<Sp_Biom(s,i-rec_age)<< " "<< " 999" << " "<< mod_rec(s,i)<<endl;
        
        R_report   << endl;

    for (r=1;r<=nreg(s);r++)
    {
      R_report <<"$SR_Curve_years_"<< (r) <<endl;
      R_report << yr_rec_est(cum_regs(s)+r) <<endl;
    }
    R_report   << endl;
    
    for (r=1;r<=nreg(s);r++)
    {
      R_report <<"$stock_Rec_Curve_"<< (r) <<endl;
      R_report <<"0 0"<<endl;
      dvariable stock;
      for (i=1;i<=300;i++)
      {
        stock = double (i) * Bzero(cum_regs(s)+r) /250.; //max(Bzero) //Bzero(cum_regs(s)+1,cum_regs(s)+nreg(s))
        if (active(log_Rzero(cum_regs(s)+r)))
          R_report << stock <<" "<< SRecruit(stock, cum_regs(s)+r)<<endl;
        else
          R_report << stock <<" 99 "<<endl;
      }
    }
    
    R_report   << endl;

    R_report   << endl<<"$Like_Comp" <<endl;
    obj_comps(13)= obj_fun - sum(obj_comps(1,12)) ; // Residual 
    obj_comps(14)= obj_fun ;                  // Total
    R_report   <<obj_comps<<endl;
    R_report   << endl;
    R_report   << endl<<"$Like_Comp_names" <<endl;
    R_report   <<"catch_like       "<<endl
             <<"age_like_fsh     "<<endl
             <<"length_like_fsh  "<<endl
             <<"sel_like_fsh     "<<endl
             <<"ind_like         "<<endl
             <<"age_like_ind     "<<endl
             <<"length_like_ind  "<<endl
             <<"sel_like_ind     "<<endl
             <<"rec_like         "<<endl
             <<"fpen             "<<endl
             <<"post_priors_indq "<<endl
             <<"post_priors      "<<endl
             <<"residual         "<<endl
             <<"total            "<<endl;
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report << "$Sel_Fshry_"<< (k) <<""<<endl;
        R_report << sel_like_fsh(k) << endl;
      }
    }
    R_report   << endl;
  
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report << "$Survey_Index_"<< (k) <<"" <<endl;
        R_report<< ind_like(k)<<endl;
      }
    }
    R_report   << endl;

    R_report << setw(10)<< setfixed() << setprecision(5) <<endl;
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report << "$Age_Survey_"<< (k) <<"" <<endl;
        R_report << age_like_ind(k)<<endl;
      }
    }
    R_report   << endl;

    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report << "$Sel_Survey_"<< (k) <<""<<endl;
        R_report<< sel_like_ind(k,1) <<" "<<sel_like_ind(k,2)<<" "<<sel_like_ind(k,3)<< endl;
      }
    }
    R_report   << endl;

    R_report << setw(10)<< setfixed() << setprecision(5) <<endl;
    R_report   << "$Rec_Pen" <<endl;
    for (r=1;r<=nreg(s);r++)
      R_report << sigmar(rec_map(s,r))<<"  ";
    R_report<<rec_like(s)<<endl;
    R_report   << endl;
    R_report<< "$m_sigmar"<<endl<<m_sigmar(cum_regs(s)+1,cum_regs(s)+nreg(s))<<endl;
    R_report<< "$sigmar"<<endl;
    for (r=1;r<=nreg(s);r++)
      R_report << sigmar(rec_map(s,r))<<"  ";
    R_report<<endl;
    R_report << endl<< "$rec_dev"<< endl;
    for (i=styr_rec;i<=endyr;i++)
      R_report << i<< " "<<rec_dev(s,i)<<endl;
    R_report   << endl;

    R_report   << "$F_Pen" <<endl;
    R_report<<fpen(1)<<"  "<<fpen(2)<<endl;
    R_report   << endl;
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report << "$Q_Survey_"<< (k) <<""<<endl
               << " "<<post_priors_indq(k)
               << " "<< q_ind(k,1)
               << " "<< qprior(k)
               << " "<< cvqprior(k)<<endl;
        R_report << "$Q_power_Survey_"<< (k) <<""<<endl
               << " "<<post_priors_indq(k)
               << " "<< q_power_ind(k)
               << " "<< q_power_prior(k)
               << " "<< cvq_power_prior(k)<<endl;
      }
    }
             R_report   << endl;
    R_report << "$Mest"<<endl;
    for (r=1;r<=nreg(s);r++)
      R_report << " "<< post_priors(mort_map(s,r),1)
               << " "<< Mest(mort_map(s,r))
               << " "<< natmortprior(mort_map(s,r))
               << " "<< cvnatmortprior(mort_map(s,r)) <<endl;
    R_report   << endl;
    R_report << "$Steep"<<endl;
    for (r=1;r<=nreg(s);r++)
      R_report << " "<< post_priors(rec_map(s,r),2)
               << " "<< steepness(rec_map(s,r))
               << " "<< steepnessprior(rec_map(s,r))
               << " "<< cvsteepnessprior(rec_map(s,r)) <<endl;
    R_report   << endl;
    R_report << "$Sigmar"<<endl;
    for (r=1;r<=nreg(s);r++)
      R_report << " "<< post_priors(rec_map(s,r),3)
               << " "<< sigmar(rec_map(s,r))
               << " "<< sigmarprior(rec_map(s,r))
               << " "<< cvsigmarprior(rec_map(s,r)) <<endl;
    R_report   << endl;
    R_report<<"$Num_parameters_Est"<<endl;
    R_report<<initial_params::nvarcalc()<<endl;
    R_report   << endl;
    
    R_report<<"$Steep_Prior" <<endl;
    for (r=1;r<=nreg(s);r++)
      R_report<<steepnessprior(rec_map(s,r))<<" "<<
        cvsteepnessprior(rec_map(s,r))<<" "<<
        phase_srec(rec_map(s,r))<<" "<< endl;
    R_report   << endl;

    R_report<<"$sigmarPrior " <<endl;
    for (r=1;r<=nreg(s);r++)
      R_report<<sigmarprior(rec_map(s,r))<<" "<<  cvsigmarprior(rec_map(s,r)) <<" "<<phase_sigmar(rec_map(s,r))<<endl;
    R_report   << endl;

    R_report<<"$Rec_estimated_in_styr_endyr " <<endl;
    R_report<<styr_rec    <<" "<<endyr        <<" "<<endl;
    R_report   << endl;
    for (r=1;r<=nreg(s);r++)
    {
      R_report<<"$SR_Curve_fit__in_styr_endyr_" << (r) <<""<<endl;
      R_report<<styr_rec_est(s,r)<<" "<<endyr_rec_est(s,r)<<" "<<endl;
    }
    R_report   << endl;
    R_report<<"$Model_styr_endyr" <<endl;
    R_report<<styr        <<" "<<endyr        <<" "<<endl;
    R_report   << endl;

    R_report<<"$M_prior "<<endl;
    for (r=1;r<=nreg(s);r++)
      R_report<< natmortprior(mort_map(s,r))<< " "<< cvnatmortprior(mort_map(s,r))<<" "<<phase_M(mort_map(s,r))<<endl;
    R_report   << endl;
    R_report<<"$qprior " <<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report<< qprior(k)<<" "<<cvqprior(k)<<" "<< phase_q(k)<<endl;
    R_report<<"$q_power_prior " <<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report<< q_power_prior(k)<<" "<<cvq_power_prior(k)<<" "<< phase_q_power(k)<<endl;
    R_report   << endl;

    R_report<<"$cv_catchbiomass " <<endl;
    R_report<<cv_catchbiomass<<" "<<endl;
    R_report   << endl;
    R_report<<"$Projection_years"<<endl;
    R_report<< nproj_yrs<<endl;
    R_report   << endl;
  
    R_report << "$Fsh_sel_opt_fish "<<endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report<<k<<" "<<fsh_sel_opt(k)<<" "<<sel_change_in_fsh(k)<<endl;
    R_report   << endl;
    R_report<<"$Survey_Sel_Opt_Survey " <<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report<<k<<" "<<(ind_sel_opt(k))<<endl;
    R_report   << endl;
    
    R_report <<"$Phase_survey_Sel_Coffs "<<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report <<phase_selcoff_ind(k)<<" "<<endl;
    R_report   << endl << endl;
    R_report <<"$Fshry_Selages " << endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report << nselages_in_fsh(k) <<" ";
    R_report   << endl << endl;
    R_report <<"$Survy_Selages " <<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report <<nselages_in_ind(k) <<" ";
    R_report   << endl << endl;

    R_report << "$Phase_for_age_spec_fishery"<<endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report <<phase_selcoff_fsh(k)<<" ";
    R_report   << endl << endl;
    R_report << "$Phase_for_logistic_fishery"<<endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report <<phase_logist_fsh(k)<<" ";
    R_report   << endl << endl;
    R_report << "$Phase_for_dble_logistic_fishery "<<endl;
    for (k=1;k<=nfsh;k++)
      if (sel_map(1,k) == s)
        R_report <<phase_dlogist_fsh(k)<<" ";
    R_report   << endl << endl;

    R_report << "$Phase_for_age_spec_survey  "<<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report <<phase_selcoff_ind(k)<<" ";
    R_report   << endl << endl;
    R_report << "$Phase_for_logistic_survey  "<<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report <<phase_logist_ind(k)<<" ";
    R_report   << endl << endl;
    R_report << "$Phase_for_dble_logistic_indy "<<endl;
    for (k=1;k<=nind;k++)
      if (sel_map(1,k+nfsh) == s)
        R_report <<phase_dlogist_ind(k)<<" ";
    R_report   << endl << endl;
  
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        if (nyrs_fsh_age(k)>0)
        {
          R_report <<"$EffN_Fsh_"<<(k)<<""<<endl;
          for (i=1;i<=nyrs_fsh_age(k);i++)
          {
            double sda_tmp = Sd_age(oac_fsh(k,i));
            R_report << yrs_fsh_age(k,i);
            R_report << " "<<Eff_N(oac_fsh(k,i),eac_fsh(k,i)) ;
            R_report << " "<<Eff_N2(oac_fsh(k,i),eac_fsh(k,i));
            R_report << " "<<mn_age(oac_fsh(k,i));
            R_report << " "<<mn_age(eac_fsh(k,i));
            R_report << " "<<sda_tmp;
            R_report << " "<<mn_age(oac_fsh(k,i)) - sda_tmp *2. / sqrt(n_sample_fsh_age(k,i));
            R_report << " "<<mn_age(oac_fsh(k,i)) + sda_tmp *2. / sqrt(n_sample_fsh_age(k,i));
            R_report <<endl;
          }
        }
      }
    }
  
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        if (nyrs_fsh_length(k)>0)
        {
          R_report <<"$EffN_Length_Fsh_"<<(k)<<""<<endl;
          for (i=1;i<=nyrs_fsh_length(k);i++)
          {
            double sda_tmp = Sd_length(olc_fsh(k,i));
            R_report << yrs_fsh_length(k,i);
            R_report << " "<<Eff_N(olc_fsh(k,i),elc_fsh(k,i)) ;
            R_report << " "<<Eff_N2_L(olc_fsh(k,i),elc_fsh(k,i));
            R_report << " "<<mn_length(olc_fsh(k,i));
            R_report << " "<<mn_length(elc_fsh(k,i));
            R_report << " "<<sda_tmp;
            R_report << " "<<mn_length(olc_fsh(k,i)) - sda_tmp *2. / sqrt(n_sample_fsh_length(k,i));
            R_report << " "<<mn_length(olc_fsh(k,i)) + sda_tmp *2. / sqrt(n_sample_fsh_length(k,i));
            R_report <<endl;
          }
        }
      }
    }


    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report <<"$C_fsh_" <<(k)<<"" << endl; 
        for (i=styr;i<=endyr;i++)
          R_report <<i<<" "<<catage(k,i)<< endl;
      }
    }

    R_report <<"$wt_a_pop" << endl<< wt_pop(s)  <<endl;
    R_report <<"$mature_a" << endl<< maturity(s)<<endl;
    for (k=1;k<=nfsh;k++)
    {
      if (sel_map(1,k) == s)
      {
        R_report <<"$wt_fsh_"<<(k)<<""<<endl;
        for (i=styr;i<=endyr;i++)
          R_report <<i<<" "<<wt_fsh(k,i)<< endl;
      }
    }
  
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        R_report <<"$wt_ind_"<<(k)<<""<<endl;
        for (i=styr;i<=endyr;i++)
          R_report <<i<<" "<<wt_ind(k,i)<< endl;
      }
    }
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        if (nyrs_ind_age(k)>0)
        {
          R_report <<"$EffN_Survey_"<<(k)<<""<<endl;
          for (i=1;i<=nyrs_ind_age(k);i++)
          {
            double sda_tmp = Sd_age(oac_ind(k,i));
            R_report << yrs_ind_age(k,i)
                     << " "<<Eff_N(oac_ind(k,i),eac_ind(k,i)) 
                     << " "<<Eff_N2(oac_ind(k,i),eac_ind(k,i))
                     << " "<<mn_age(oac_ind(k,i))
                     << " "<<mn_age(eac_ind(k,i))
                     << " "<<sda_tmp
                     << " "<<mn_age(oac_ind(k,i)) - sda_tmp *2. / sqrt(n_sample_ind_age(k,i))
                     << " "<<mn_age(oac_ind(k,i)) + sda_tmp *2. / sqrt(n_sample_ind_age(k,i))
                     <<endl;
          }
        }
      }
    }
    for (k=1;k<=nind;k++)
    {
      if (sel_map(1,k+nfsh) == s)
      {
        if (nyrs_ind_length(k)>0)
        {
          R_report <<"$EffN_Length_Survey_"<<(k)<<""<<endl;
          for (i=1;i<=nyrs_ind_length(k);i++)
          {
            double sda_tmp = Sd_age(olc_ind(k,i));
            R_report << yrs_ind_length(k,i)
                     << " "<<Eff_N(olc_ind(k,i),elc_ind(k,i)) 
                     << " "<<Eff_N2_L(olc_ind(k,i),elc_ind(k,i))
                     << " "<<mn_length(olc_ind(k,i))
                     << " "<<mn_length(elc_ind(k,i))
                     << " "<<sda_tmp
                     << " "<<mn_length(olc_ind(k,i)) - sda_tmp *2. / sqrt(n_sample_ind_length(k,i))
                     << " "<<mn_length(olc_ind(k,i)) + sda_tmp *2. / sqrt(n_sample_ind_length(k,i))
                     <<endl;
          }
        }
      }
    }
  
    R_report<<"$msy_mt"<<endl; 
    dvar_matrix sel_tmp(1,nages,1,nfsh);
    dvar_vector sumF(1,nstk);
    sel_tmp.initialize();
    for (i=styr;i<=endyr;i++) 
    { 
      sumF.initialize();
      for (k=1;k<=nfsh;k++)
      {
        Fratio(k) = sum(F(k,i)) ;
        sumF(sel_map(1,k)) += Fratio(k) ;
      }
      for (k=1;k<=nfsh;k++)
        Fratio(k) /= sumF(sel_map(1,k));
      sumF /= nages;
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == s)
          for (j=1;j<=nages;j++)
            sel_tmp(j,k) = sel_fsh(k,i,j); 
      get_msy(i);
      // important for time-varying natural mortality...
      dvariable spr_mt_ft = spr_ratio(sumF(s),sel_tmp,i,s)  ;
      // Yr Fspr 1-Fspr F/Fmsy Fmsy F Fsprmsy MSY MSYL Bmsy Bzero SSB B/Bmsy
      R_report<< i<<
              " "<< spr_mt_ft                     <<
              " "<< (1.-spr_mt_ft)                << 
              " "<< Fcur_Fmsy(s)                  <<
              " "<< Fmsy(s)                       <<
              " "<< sumF(s)                       <<
              " "<< spr_ratio(Fmsy(s),sel_tmp,i,s)<<
              " "<< MSY(s)                        <<
              " "<< MSYL(s)                       <<
              " "<< Bmsy(s)                       <<
              " "<< Bzero(cum_regs(s)+yy_sr(s,i)) <<
              " "<< Sp_Biom(s,i)                  <<
              " "<< Bcur_Bmsy(s)                  <<
              endl ;
    }
    for (r=1;r<=nreg(s);r++)
    {
      R_report<<"$age2len_"<<growth_map(s,r)<<endl;
      for (j=1;j<=nages;j++)
        R_report<<P_age2len(growth_map(s,r),j)<<endl;
    }
    R_report<<"$msy_m0"<<endl; 
    sel_tmp.initialize();
    // NOTE Danger here
    dvar3_array mtmp = M;
    for (i=styr;i<=endyr;i++) 
    { 
      M(s,i) = M(s,styr);
      sumF.initialize();
      for (k=1;k<=nfsh;k++)
      {
        Fratio(k) = sum(F(k,i)) ;
        sumF(sel_map(1,k)) += Fratio(k) ;
      }
      for (k=1;k<=nfsh;k++)
        Fratio(k) /= sumF(sel_map(1,k));
      for (k=1;k<=nfsh;k++)
        if (sel_map(1,k) == s)
          for (j=1;j<=nages;j++)
            sel_tmp(j,k) = sel_fsh(k,i,j); 
      get_msy(i);
      sumF /= nages;
      // important for time-varying natural mortality...
      dvariable spr_mt_ft = spr_ratio(sumF(s),sel_tmp,i,s)  ;
      dvariable spr_mt_f0 = spr_ratio(0.,sel_tmp,i,s)  ;
      R_report<< i<<
              " "<< spr_mt_ft                     <<
              " "<< spr_mt_f0                     <<
              " "<< (1.-spr_mt_f0)/(1-spr_mt_ft)  << 
              " "<< Fcur_Fmsy(s)                  <<
              " "<< Fmsy(s)                       <<
              " "<< sumF(s)                       <<
              " "<< spr_ratio(Fmsy(s),sel_tmp,i,s)<<
              " "<< MSY(s)                        <<
              " "<< Bmsy(s)                       <<
              " "<< MSYL(s)                       <<
              " "<< Bcur_Bmsy(s)                  <<
              endl ;
    }
    // M = mtmp;
    R_report<< "$F40_est"<<endl<<F40_est(s)<<endl;
    R_report<< "$F35_est"<<endl<<F35_est(s)<<endl;

		// Save components of recruitment likelihood to R
    for (k=1;k<=nstk;k++)
		{
      R_report   << "$rec_like_SRR_stock_" <<k<<endl<<rec_like(k,1)<< endl;
      R_report   << "$rec_like_pen_stock_" <<k<<endl<<rec_like(k,2)<< endl;
      R_report   << "$rec_like_fut_stock_" <<k<<endl<<rec_like(k,3)<< endl;
		}
    R_report.close();
  }


FUNCTION double mn_age(const dvector& pobs)
  // int lb1 = pobs.indexmin();
  // int ub1 = pobs.indexmax();
  // dvector av = age_vector(lb1,ub1)  ;
  // double mobs = value(pobs.shift(rec_age)*age_vector);
  double mobs = (pobs*age_vector);
  return mobs;

FUNCTION double mn_age(const dvar_vector& pobs)
  // int lb1 = pobs.indexmin();
  // int ub1 = pobs.indexmax();
  // dvector av = age_vector(lb1,ub1)  ;
  // double mobs = value(pobs.shift(rec_age)*age_vector);
  double mobs = value(pobs*age_vector);
  return mobs;

FUNCTION double Sd_age(const dvector& pobs)
  // double mobs = (pobs.shift(rec_age)*age_vector);
  // double stmp = (sqrt(elem_prod(age_vector,age_vector)*pobs.shift(rec_age) - mobs*mobs));
  double mobs = (pobs*age_vector);
  double stmp = sqrt((elem_prod(age_vector,age_vector)*pobs) - mobs*mobs);
  return stmp;

FUNCTION double mn_length(const dvector& pobs)
  double mobs = (pobs*len_bins);
  return mobs;

FUNCTION double mn_length(const dvar_vector& pobs)
  double mobs = value(pobs*len_bins);
  return mobs;

FUNCTION double Sd_length(const dvector& pobs)
  double mobs = (pobs*len_bins);
  double stmp = sqrt((elem_prod(len_bins,len_bins)*pobs) - mobs*mobs);
  return stmp;

FUNCTION double Eff_N_adj(const double, const dvar_vector& pobs, const dvar_vector& phat)
  int lb1 = pobs.indexmin();
  int ub1 = pobs.indexmax();
  dvector av = age_vector(lb1,ub1)  ;
  double mobs = value(pobs*av);
  double mhat = value(phat*av );
  double rtmp = mobs-mhat;
  double stmp = value(sqrt(elem_prod(av,av)*pobs - mobs*mobs));
  return square(stmp)/square(rtmp);

FUNCTION double Eff_N2(const dvector& pobs, const dvar_vector& phat)
  int lb1 = pobs.indexmin();
  int ub1 = pobs.indexmax();
  dvector av = age_vector(lb1,ub1)  ;
  double mobs =      (pobs*av);
  double mhat = value(phat*av );
  double rtmp = mobs-mhat;
  double stmp = (sqrt(elem_prod(av,av)*pobs - mobs*mobs));
  return square(stmp)/square(rtmp);

FUNCTION double Eff_N(const dvector& pobs, const dvar_vector& phat)
  dvar_vector rtmp = elem_div((pobs-phat),sqrt(elem_prod(phat,(1-phat))));
  double vtmp;
  vtmp = value(norm2(rtmp)/size_count(rtmp));
  return 1./vtmp;

FUNCTION double Eff_N2_L(const dvector& pobs, const dvar_vector& phat)
  dvector av = len_bins  ;
  double mobs =      (pobs*av);
  double mhat = value(phat*av );
  double rtmp = mobs-mhat;
  double stmp = (sqrt(elem_prod(av,av)*pobs - mobs*mobs));
  return square(stmp)/square(rtmp);

FUNCTION double get_AC(const int& indind)
  // Functions to compute autocorrelation in residuals 
  int i1,i2,iyr;
  i1 = 1;
  i2 = nyrs_ind(indind);
  double actmp;
  dvector res(1,i2);
  for (i=1;i<=i2;i++)
  {
    iyr = int(yrs_ind(indind,i));
    res(i) = log(obs_ind(indind,i)) - value(log(pred_ind(indind,i)));
  }
  double m1 = (mean(res(i1,i2-1)));
  double m2 = (mean(res(i1+1,i2))); 
  actmp = mean( elem_prod( ++res(i1,i2-1) - m1, res(i1+1,i2) - m2)) /
          (sqrt(mean( square(res(i1,i2-1) - m1 )))  * sqrt(mean(square(res(i1+1,i2) - m2 ))) );
  return(actmp);

FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred)+0.000001;
  sdnr = std_dev(elem_div(obs+0.000001-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;

  /**
   * @brief Calculate Francis weights
   * @details this code based on equation TA1.8 in Francis(2011) should be changed so separate weights if by sex
   *
   * Produces the new weight that should be used.
  **/
FUNCTION double calc_Francis_weights(const dmatrix oac, const dvar_matrix eac, const ivector sam )
  {
    int nobs;
    int i1=oac.rowmin();
    int i2=oac.rowmax();
    double lfwt,Var,Pre,Obs;
    dvector ages(oac.colmin(),nages);
    for (int i=oac.colmin();i<=nages;i++) 
      ages(i) = double(i)+.5;
    nobs = oac.rowsize();
    dvector resid(i1,i2);
    resid.initialize();
    for ( int i = i1; i <= i2; i++ )
    {
      // Obs = sum(elem_prod(oac(i), ages+.5));
      Obs = oac(i) * (ages+.5);
      // Pre = sum(elem_prod(value(eac(i)), ages+.5));
      Pre = value(eac(i)) * (ages+.5);
      Var = value(eac(i)) * square(ages+.5);
      Var -= square(Pre);
      resid(i) = (Obs - Pre) / sqrt(Var * 1.0 / (sam(i) ));
    }
    lfwt = 1.0 / (square(std_dev(resid)) * ((nobs - 1.0) / nobs * 1.0));
    // lfwt(k) *= lf_lambda(k);
    return lfwt;
  }



