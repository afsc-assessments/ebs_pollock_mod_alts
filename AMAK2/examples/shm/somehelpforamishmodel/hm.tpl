///////////////////////////////////////////////////////////////////
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
//    obs_se_ind    Observed standard errors (annual)
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
//
//////////////////////////////////////////////////////////////////////////////
 // To ADD/FIX:
 //   parameterization of steepness to work the same (wrt prior) for ricker and bholt
 //   natural mortality-at-age
 //   splines for selectivity
//////////////////////////////////////////////////////////////////////////////

DATA_SECTION
  !!version_info+="HorseMackerel_2011";
  int iseed 
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
  writeinput<<version_info<<endl;
  tmpstring=adprogram_name + adstring(".dat");
  // if (argc > 1)
  // {
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
      cout<<"Got to operating model option "<<oper_mod<<endl;
    }
    if ( (on=option_match(argc,argv,"-mcmc"))>-1)
    {
      mcmcmode = 1;
    }
  // }
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
  int nages
  !!  nages = oldest_age - rec_age + 1;
  int styr_rec
  int styr_sp
  int endyr_sp
  int nyrs
  !!  nyrs = endyr - styr + 1;
  !! styr_rec = (styr - nages) + 1;    // First year of recruitment
  !! styr_sp  = styr_rec - rec_age ;   // First year of spawning biomass  
  !! endyr_sp = endyr   - rec_age ;    // endyr year of (main) spawning biomass
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
  init_matrix catch_bio(1,nfsh,styr,endyr)
  init_matrix catch_bio_sd(1,nfsh,styr,endyr)   // Specify catch-estimation precision
  !! log_input(catch_bio);
  !! log_input(catch_bio_sd);
  matrix catch_bio_lsd(1,nfsh,styr,endyr)      //Catch biomass standard errors (for lognormal)
  matrix catch_bio_lva(1,nfsh,styr,endyr)      //Catch biomass standard errors (for lognormal)
  !! catch_bio_lsd = sqrt(log(square(catch_bio_sd) + 1.));
  !! catch_bio_lva = log(square(catch_bio_sd) + 1.);
  matrix catch_bioT(styr,endyr,1,nfsh)
  !! catch_bioT = trans(catch_bio);
  // !! cout << catch_bio<<endl;exit(1);

  init_ivector nyrs_fsh_age(1,nfsh)
  !! log_input(nyrs_fsh_age);

  init_imatrix yrs_fsh_age(1,nfsh,1,nyrs_fsh_age)
  !! log_input(yrs_fsh_age);

  init_matrix n_sample_fsh_age(1,nfsh,1,nyrs_fsh_age)    //Years of index value (annual)
  !! log_input(n_sample_fsh_age);

  init_3darray oac_fsh(1,nfsh,1,nyrs_fsh_age,1,nages)
    3darray sim_catage(1,nfsh,1,nyrs_fsh_age,1,nages);
  !! log_input(oac_fsh);

  init_3darray wt_fsh(1,nfsh,styr,endyr,1,nages)  //values of index proportions at age
  !! log_input(wt_fsh);

//  Define indices
  init_int nind                                   //number of indices
  !! log_input(nind);
  int nfsh_and_ind
  !! nfsh_and_ind = nfsh+nind;
  imatrix pindname(1,nind,1,2)
  init_adstring indnameread;
 LOCAL_CALCS
  for(k=1;k<=nind;k++) 
  {
    pindname(k,1)=1; 
    pindname(k,2)=1;
  }    // set whole array to equal 1 in case not enough names are read
  k=1;
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
  !! log_input(nyrs_ind);
  init_imatrix yrs_ind(1,nind,1,nyrs_ind)         //Years of index value (annual)
  !! log_input(yrs_ind);
  init_vector mo_ind(1,nind)                      //Month indices occur 
  !! log_input(mo_ind);
  vector ind_month_frac(1,nind)
  !! ind_month_frac = (mo_ind-1.)/12.;
  init_matrix obs_ind(1,nind,1,nyrs_ind)          //values of indices index value (annual)
  !! log_input(obs_ind);
  init_matrix obs_se_ind(1,nind,1,nyrs_ind)       //values of indices serrs
  !! log_input(obs_se_ind);
  matrix     obs_lse_ind(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
  matrix     obs_lva_ind(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
  matrix        corr_dev(1,nind,1,nyrs_ind) //Index standard errors (for lognormal)
  matrix        corr_eff(1,nfsh,styr,endyr) //Index standard errors (for lognormal)
  matrix         act_eff(1,nfsh,styr,endyr) //Index standard errors (for lognormal)
  vector              ac(1,nind);
  !! obs_lse_ind = elem_div(obs_se_ind,obs_ind);
  !! obs_lse_ind = sqrt(log(square(obs_lse_ind) + 1.));
  !! obs_lva_ind = square(obs_lse_ind);

  init_ivector nyrs_ind_age(1,nind)               //Number of years of indices index value (annual)
  !! log_input(nyrs_ind_age);
  init_matrix yrs_ind_age(1,nind,1,nyrs_ind_age)  //Years of indices index value (annual)
  !! log_input(yrs_ind_age);
  init_matrix n_sample_ind_age(1,nind,1,nyrs_ind_age)         //Years of indices index value (annual)
  !! log_input(n_sample_ind_age);
  init_3darray oac_ind(1,nind,1,nyrs_ind_age,1,nages)  //values of Index proportions at age
  !! log_input(oac_ind);
  init_3darray  wt_ind(1,nind,styr,endyr,1,nages)      //values of Index proportions at age
  !! log_input(wt_ind);

  vector age_vector(1,nages);
    !! for (int j=1;j<=nages;j++)
      !!  age_vector(j) = double(j+rec_age-1);
  init_vector wt_pop(1,nages)
  !! log_input(wt_pop);
  init_vector maturity(1,nages)
   !! if (max(maturity)>.9) maturity /=2.;
   !! log_input(maturity);
  vector wt_mature(1,nages);
  !! wt_mature = elem_prod(wt_pop,maturity) ;

  //Spawning month-----
  init_number spawnmo
  number spmo_frac
  !! spmo_frac = (spawnmo-1)/12.;
  !! log_input(spawnmo);
  !! log_input(spmo_frac);

  init_matrix age_err(1,nages,1,nages)
  !! log_input(age_err);

  int k // Index for fishery or indices
  int i // Index for year
  int j // Index for age
 LOCAL_CALCS
  // Rename data file to the control data section... 
  ad_comm::change_datafile_name(cntrlfile_name);
  *(ad_comm::global_datafile) >>  datafile_name; 
  *(ad_comm::global_datafile) >>  model_name; 
  log_input(cntrlfile_name);
 END_CALCS
  // Matrix of selectivity mappings--row 1 is type (1=fishery, 2=index) and row 2 is index within that type
  //  e.g., the following for 2 fisheries and 4 indices means that index 3 uses fishery 1 selectivities,
  //         the other fisheries and indices use their own parameterization
  //  1 1 2 2 1 2 
  //  1 2 1 2 1 4
  init_imatrix sel_map(1,2,1,nfsh_and_ind) 
  // maps fisheries and indices into sequential sel_map for sharing purposes
  !! writeinput<< "# Map shared selectivity: "<< endl;log_input(sel_map);
  !! log_input(datafile_name);
  !! log_input(model_name);
  !! projfile_name = cntrlfile_name(1,length(cntrlfile_name)-4) + ".prj";
  init_int    SrType        // 2 Bholt, 1 Ricker
  init_int use_age_err      // nonzero value means use...
  init_int use_effort       // Effort time series (not implemented fully)
  init_number steepnessprior
  init_number cvsteepnessprior
  init_int    phase_srec

  init_number sigmarprior
  number log_sigmarprior
  init_number cvsigmarprior
  init_int    phase_sigmar
  init_int    styr_rec_est
  init_int    endyr_rec_est
  int nrecs_est;

  init_number natmortprior
  init_number cvnatmortprior
  init_int    phase_M
  init_vector M_offset(1,nages) 

  init_vector qprior(1,nind)      
  vector log_qprior(1,nind)      
  init_vector cvqprior(1,nind)     
  init_ivector phase_q(1,nind)

  init_vector q_power_prior(1,nind)      
  vector log_q_power_prior(1,nind)      
  init_vector cvq_power_prior(1,nind)     
  init_ivector phase_q_power(1,nind)

  init_ivector    q_age_min(1,nind)     // Age that q relates to...
  init_ivector    q_age_max(1,nind)     // Age that q relates to...

  init_number cv_catchbiomass
  number catchbiomass_pen
  !!catchbiomass_pen= 1./(2*cv_catchbiomass*cv_catchbiomass);
  init_int nproj_yrs

  int styr_fut
  int endyr_fut            // LAst year for projections
  int phase_Rzero
  int phase_nosr
  number Steepness_UB
  !! phase_Rzero =  4;
  !! phase_nosr  = -3;
  !! styr_fut    = endyr+1;
  !! endyr_fut   = endyr + nproj_yrs; 

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
  matrix sel_change_in_fsh(1,nfsh,styr,endyr);
  ivector nselages_in_fsh(1,nfsh)

  ivector n_sel_ch_fsh(1,nfsh);
  ivector n_sel_ch_ind(1,nind);
  imatrix yrs_sel_ch_tmp(1,20,1,endyr-styr+1);
  imatrix yrs_sel_ch_tmp_ind(1,nind,1,endyr-styr+1);
  !! yrs_sel_ch_tmp.initialize();
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
  sel_change_in_fsh.initialize()   ;  
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
  for (k=1;k<=nfsh;k++)
  {
    *(ad_comm::global_datafile) >> fsh_sel_opt(k)  ;  
    writeinput<<fshname(k)<<" fishery Sel Option: "<<fsh_sel_opt(k)<<endl;
    switch (fsh_sel_opt(k))
    {
      case 1 : // Selectivity coefficients 
      {
        *(ad_comm::global_datafile) >> nselages_in_fsh(k)   ;  
        *(ad_comm::global_datafile) >> phase_sel_fsh(k);  
        *(ad_comm::global_datafile) >> curv_pen_fsh(k) ;
        *(ad_comm::global_datafile) >> seldec_pen_fsh(k) ;
        seldec_pen_fsh(k) *= seldec_pen_fsh(k) ;
        for (int i=styr;i<=endyr;i++)
          *(ad_comm::global_datafile) >> sel_change_in_fsh(k,i) ;

        sel_change_in_fsh(k,styr)=1.; 
        n_sel_ch_fsh(k)=0.;
        int j=1; 
        yrs_sel_ch_tmp(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_fsh(k,i)>0) {
            j++; yrs_sel_ch_tmp(k,j) = i; } }
        n_sel_ch_fsh(k) = j; 

        // This to read in pre-specified selectivity values...
        sel_fsh_tmp.initialize();
        log_selcoffs_fsh_in.initialize();
        for (j=1;j<=nages;j++) 
            *(ad_comm::global_datafile) >> sel_fsh_tmp(j);  
        for (int jj=2;jj<=n_sel_ch_fsh(k);jj++) 
        {
          // Set the selectivity for the oldest group
          for (j=nselages_in_fsh(k)+1;j<=nages;j++) 
          {
            sel_fsh_tmp(j)  = sel_fsh_tmp(nselages_in_fsh(k));  
          }
          // Set tmp to actual initial vectors...
          log_selcoffs_fsh_in(k,jj)(1,nselages_in_fsh(k)) = log(sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7/mean(sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7) );
          writeinput<<"Sel_in_fsh "<< mfexp(log_selcoffs_fsh_in(k,jj))<<endl;
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
        for (int i=styr;i<=endyr;i++) { 
          *(ad_comm::global_datafile) >> sel_change_in_fsh(k,i) ; }

        sel_change_in_fsh(k,styr)=1.; 
        n_sel_ch_fsh(k)=0;
        int j=1; 
        yrs_sel_ch_tmp(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_fsh(k,i)>0) {
            j++; yrs_sel_ch_tmp(k,j) = i; } }
        n_sel_ch_fsh(k) = j; 

        // This to read in pre-specified selectivity values...
        *(ad_comm::global_datafile) >> sel_slp_in_fsh(k,1) ;
        *(ad_comm::global_datafile) >> sel_inf_in_fsh(k,1) ;
        for (int jj=2;jj<=n_sel_ch_fsh(k);jj++) 
        {
          sel_inf_in_fsh(k,jj)    =     sel_inf_in_fsh(k,1) ;
          logsel_slp_in_fsh(k,jj) = log(sel_slp_in_fsh(k,1)) ;
        }

        phase_selcoff_fsh(k) = -1;
        phase_logist_fsh(k) = phase_sel_fsh(k);
        phase_dlogist_fsh(k) = -1;

        logsel_slp_in_fshv(k) = logsel_slp_in_fsh(k,1);
           sel_inf_in_fshv(k) =    sel_inf_in_fsh(k,1);
        break;
      }
      case 3 : // Double logistic 
      {
        writeinput << "Double logistic abandoned..."<<endl;exit(1);
        *(ad_comm::global_datafile) >> phase_sel_fsh(k);  
        for (int i=styr;i<=endyr;i++) { 
          *(ad_comm::global_datafile) >> sel_change_in_fsh(k,i) ; }
        sel_change_in_fsh(k,styr)=1.; 
        n_sel_ch_fsh(k)=0;
        int j=1; 
        yrs_sel_ch_tmp(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_fsh(k,i)>0) {
            j++; yrs_sel_ch_tmp(k,j) = i; } }
        n_sel_ch_fsh(k) = j; 

        // This to read in pre-specified selectivity values...
        for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
        {
          *(ad_comm::global_datafile) >> sel_slp_in_fsh(k,jj) ;
          *(ad_comm::global_datafile) >> sel_inf_in_fsh(k,jj) ;
          *(ad_comm::global_datafile) >> sel_dslp_in_fsh(k,jj) ;
          *(ad_comm::global_datafile) >> sel_dinf_in_fsh(k,jj) ;
          logsel_slp_in_fsh(k,jj) = log(sel_slp_in_fsh(k,jj)) ;
          logsel_dslp_in_fsh(k,jj) = log(sel_dslp_in_fsh(k,jj)) ;
        }


        phase_selcoff_fsh(k) = -1;
        phase_logist_fsh(k)  = phase_sel_fsh(k);
        phase_dlogist_fsh(k) = phase_sel_fsh(k)+1;


        logsel_slp_in_fshv(k) = logsel_slp_in_fsh(k,1);
           sel_inf_in_fshv(k) =    sel_inf_in_fsh(k,1);
        logsel_dslp_in_fshv(k) = logsel_dslp_in_fsh(k,1);
           sel_dinf_in_fshv(k) =    sel_dinf_in_fsh(k,1);
        break;
      }
      case 4 : // Splines         
      {
        /*
        writeinput << "Spline selectivity option for fishery "<<fshname(k)<<endl;
        *(ad_comm::global_datafile) >> nselages_in_fsh(k)   ; // Number of selected ages 
        *(ad_comm::global_datafile) >> phase_sel_fsh(k);       
        writeinput << "  number of sel ages: "<<nselages_in_fsh(k)<<endl;
        *(ad_comm::global_datafile) >> nnodes_fsh(k)   ;      // Number of nodes, to run over selected ages only
        writeinput << "  number of nodes : "<<nnodes_fsh(k)<<endl;
        // Read in year-year changes in selectivity
        for (int i=styr;i<=endyr;i++)
          *(ad_comm::global_datafile) >> sel_change_in_fsh(k,i) ;

        sel_change_in_fsh(k,styr)=1.; n_sel_ch_fsh(k)=0;
        int j=1; 
        yrs_sel_ch_tmp(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_fsh(k,i)>0) {
            j++; yrs_sel_ch_tmp(k,j) = i; } }
        n_sel_ch_fsh(k) = j; 

        // This to read in pre-specified selectivity values...
        sel_fsh_tmp.initialize();
        log_sel_spl_fsh_in.initialize();
        for (j=1;j<=nnodes_fsh(k);j++) 
            *(ad_comm::global_datafile) >> sel_fsh_tmp(j);  

        log_sel_spl_fsh_in(k,1)(1,nnodes_fsh(k)) = log(sel_fsh_tmp(1,nnodes_fsh(k)));
        writeinput<<"initial spline nodal values "<< mfexp(log_sel_spl_fsh_in(k,1))<<endl;
        // exit(1);
        phase_sel_spl_fsh(k) = phase_sel_fsh(k);
        phase_selcoff_fsh(k) = -1;
        phase_logist_fsh(k)  = -1;
        phase_dlogist_fsh(k) = -1;
        //  xnodes increments from 0-1 by number of nodes
        xnodes_fsh(k).fill_seqadd(0,1.0/(nnodes_fsh(k)-1));
        // xages_fsh increments from 0-1 by number of ages, say
        xages_fsh(k).fill_seqadd(0,1.0/(nages-1));
        // xages_fsh(k).fill_seqadd(0,1.0/(nselages_in_fsh(k)-1)); //prefer to use nselages but need 3d version to work
        */
      }
      break;
      writeinput << fshname(k)<<" fish sel opt "<<endl<<fsh_sel_opt(k)<<" "<<endl<<"Sel_change"<<endl<<sel_change_in_fsh(k)<<endl;
    }
  }
  // Index here..............
  for(k=1;k<=nind;k++)
  {
    *(ad_comm::global_datafile) >> ind_sel_opt(k)  ;  
    writeinput << endl<<"Index "<<indname(k)<<" Selectivity option: "<<ind_sel_opt(k)  <<endl;  
    switch (ind_sel_opt(k))
    {
      case 1 : // Selectivity coefficients  indices
      {
        *(ad_comm::global_datafile) >> nselages_in_ind(k)   ;  
        *(ad_comm::global_datafile) >> phase_sel_ind(k);  
        *(ad_comm::global_datafile) >> curv_pen_ind(k) ;
        *(ad_comm::global_datafile) >> seldec_pen_ind(k) ;
        seldec_pen_ind(k) *= seldec_pen_ind(k) ;
        writeinput << "Index: "          <<indname(k)<<endl; // cin>>junk;
        writeinput << "Index selages: "  <<nselages_in_ind(k)<<endl; // cin>>junk;
        writeinput << "Index sel_phase: "<< phase_sel_ind(k)<<endl; // cin>>junk;
        writeinput << "ind seldec_penalt: "<<seldec_pen_ind(k) << endl;
        for (int i=styr;i<=endyr;i++)
          *(ad_comm::global_datafile) >> sel_change_in_ind(k,i) ;

        sel_change_in_ind(k,styr)=1.; 
        n_sel_ch_ind(k)=0; 
        int j=1; 
        yrs_sel_ch_tmp_ind(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_ind(k,i)>0) {
            j++; yrs_sel_ch_tmp_ind(k,j) = i; } }
        n_sel_ch_ind(k) = j;
        writeinput << "ind sel_change_in: "<<sel_change_in_ind(k) << endl;
        writeinput <<"Number of select changes: indices "<<k<<": "<<n_sel_ch_ind(k)<<endl;
        writeinput<<"Year sel change 1 "<<yrs_sel_ch_tmp_ind(1,n_sel_ch_ind(k))<<endl;// exit(1);

        // This to read in pre-specified selectivity values...
        sel_ind_tmp.initialize();
        log_selcoffs_ind_in.initialize();
        for (j=1;j<=nages;j++) 
            *(ad_comm::global_datafile) >> sel_ind_tmp(j);  
        for (int jj=2;jj<=n_sel_ch_ind(k);jj++) 
        {
          for (j=nselages_in_ind(k)+1;j<=nages;j++) 
          {
            sel_ind_tmp(j)  = sel_ind_tmp(nselages_in_ind(k));  
          }
          // Set tmp to actual initial vectors...
          log_selcoffs_ind_in(k,jj)(1,nselages_in_ind(k)) = log(sel_ind_tmp(1,nselages_in_ind(k))+1e-7/mean(sel_fsh_tmp(1,nselages_in_fsh(k))+1e-7) );
          writeinput<<"Sel_in_ind "<< mfexp(log_selcoffs_ind_in(k,jj))<<endl;
        }
        phase_selcoff_ind(k) = phase_sel_ind(k);
        phase_logist_ind(k)  = -1;
        phase_dlogist_ind(k)  = -1;
      }

      break;

      case 2 : // Single logistic
      {
        *(ad_comm::global_datafile) >> phase_sel_ind(k);  
        // cout<<"PHase sel ind "<<phase_sel_ind(k)<<endl;

        for (int i=styr;i<=endyr;i++) { 
          *(ad_comm::global_datafile) >> sel_change_in_ind(k,i) ; }

        sel_change_in_ind(k,styr)=1.; n_sel_ch_ind(k)=0; int j=1; 
        writeinput << "ind sel_change_in: "<<sel_change_in_ind(k) << endl;
        yrs_sel_ch_tmp_ind(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) { if(sel_change_in_ind(k,i)>0) { j++; yrs_sel_ch_tmp_ind(k,j) = i; } }
        n_sel_ch_ind(k) = j; writeinput <<"Number of select changes: index "<<k<<": "<<n_sel_ch_ind(k)<<endl;
        writeinput<<"Year sel change 1 "<<yrs_sel_ch_tmp_ind(1,n_sel_ch_ind(k))<<endl;// exit(1);

        // This to read in pre-specified selectivity values...
        *(ad_comm::global_datafile) >> sel_slp_in_ind(k,1) ;
        *(ad_comm::global_datafile) >> sel_inf_in_ind(k,1) ;
        for (int jj=2;jj<=n_sel_ch_ind(k);jj++) {
          sel_inf_in_ind(k,jj)    =     sel_inf_in_ind(k,1) ;
          logsel_slp_in_ind(k,jj) = log(sel_slp_in_ind(k,1)) ;
        }

        writeinput<<k<<" "<< phase_sel_ind(k)<<endl;
        writeinput<<k<<" "<< sel_slp_in_ind(k)(1,n_sel_ch_ind(k)) <<endl;// exit(1);
        phase_selcoff_ind(k) = -1;
        phase_logist_ind(k) = phase_sel_ind(k);
        phase_dlogist_ind(k)  = -1;

        logsel_slp_in_indv(k) = logsel_slp_in_ind(k,1);
           sel_inf_in_indv(k) =    sel_inf_in_ind(k,1);
      }
      break;
      case 3 : // Double logistic 
      {
        writeinput << "Double logistic abandoned..."<<endl;exit(1);
        *(ad_comm::global_datafile) >> phase_sel_ind(k);  
        for (int i=styr;i<=endyr;i++) { 
          *(ad_comm::global_datafile) >> sel_change_in_ind(k,i) ; }
        sel_change_in_ind(k,styr)=1.; n_sel_ch_ind(k)=0; int j=1; 
        writeinput << "ind sel_change_in: "<<sel_change_in_ind(k) << endl;
        yrs_sel_ch_tmp_ind(k,j) = styr;
       // Number of selectivity changes is equal to the number of vectors (yr 1 is baseline)
        for (int i=styr+1;i<=endyr;i++) {
          if(sel_change_in_ind(k,i)>0) {
            j++; yrs_sel_ch_tmp_ind(k,j) = i; } }
        n_sel_ch_ind(k) = j; writeinput <<"Number of select changes: index "<<k<<": "<<n_sel_ch_ind(k)<<endl;
        writeinput<<"Year sel change 1 "<<yrs_sel_ch_tmp_ind(1,n_sel_ch_ind(k))<<endl;// exit(1);

        // This to read in pre-specified selectivity values...
        for (int jj=1;jj<=n_sel_ch_ind(k);jj++) {
          *(ad_comm::global_datafile) >> sel_slp_in_ind(k,jj) ;
          *(ad_comm::global_datafile) >> sel_inf_in_ind(k,jj) ;
          *(ad_comm::global_datafile) >> sel_dslp_in_ind(k,jj) ;
          *(ad_comm::global_datafile) >> sel_dinf_in_ind(k,jj) ;
          logsel_slp_in_ind(k,jj) = log(sel_slp_in_ind(k,jj)) ;
          logsel_dslp_in_ind(k,jj) = log(sel_dslp_in_ind(k,jj)) ;
        }


        phase_selcoff_ind(k) = -1;
        phase_logist_ind(k)  = phase_sel_ind(k);
        phase_dlogist_ind(k) = phase_sel_ind(k)+1;

        logsel_slp_in_indv(k) = logsel_slp_in_ind(k,1);
           sel_inf_in_indv(k) =    sel_inf_in_ind(k,1);
        logsel_dslp_in_indv(k) = logsel_dslp_in_ind(k,1);
           sel_dinf_in_indv(k) =    sel_dinf_in_ind(k,1);
      }
        break;
      case 4 : // spline for index
      {
      }
      break;
    }
    writeinput << indname(k)<<" ind sel opt "<<ind_sel_opt(k)<<" "<<sel_change_in_ind(k)<<endl;
  }
  writeinput<<"Phase index Sel_Coffs: "<<phase_selcoff_ind<<endl; 
 END_CALCS
  init_number test;
  !! writeinput<<" Test: "<<test<<endl;
 !! if (test!=123456789) {cerr<<"Control file not read in correctly... "<<endl;exit(1);}


  ivector nopt_fsh(1,2) // number of options...
  !! nopt_fsh.initialize();
  !! for (k=1;k<=nfsh;k++) if(fsh_sel_opt(k)==1) nopt_fsh(1)++;else nopt_fsh(2)++;

  // Fishery selectivity description:
  // type 1
  
  // Number of ages

  !! writeinput << "Fshry Selages: " << nselages_in_fsh  <<endl;
  !! writeinput << "Srvy  Selages: " << nselages_in_ind <<endl;



  !! writeinput << "Phase for age-spec fishery "<<phase_selcoff_fsh<<endl;
  !! writeinput << "Phase for logistic fishery "<<phase_logist_fsh<<endl;
  !! writeinput << "Phase for dble logistic fishery "<<phase_dlogist_fsh<<endl;

  !! writeinput << "Phase for age-spec index  "<<phase_selcoff_ind<<endl;
  !! writeinput << "Phase for logistic index  "<<phase_logist_ind<<endl;
  !! writeinput << "Phase for dble logistic indy "<<phase_dlogist_ind<<endl;

  !! for (k=1;k<=nfsh;k++) if (phase_selcoff_fsh(k)>0) curv_pen_fsh(k) = 1./ (square(curv_pen_fsh(k))*2);
  !! writeinput<<"Curv_pen_fsh: "<<curv_pen_fsh<<endl;
  !! for (k=1;k<=nind;k++) if (phase_selcoff_ind(k)>0) curv_pen_ind(k) = 1./ (square(curv_pen_ind(k))*2);

  int  phase_fmort;
  int  phase_proj;
  imatrix yrs_sel_ch_fsh(1,nfsh,1,n_sel_ch_fsh);
  imatrix   nselages_fsh(1,nfsh,1,n_sel_ch_fsh);
  matrix xnodes_fsh(1,nfsh,1,nnodes_fsh)
  matrix xages_fsh(1,nfsh,1,nages)

  imatrix yrs_sel_ch_ind(1,nind,1,n_sel_ch_ind);
  imatrix   nselages_ind(1,nind,1,n_sel_ch_ind);
 LOCAL_CALCS
  for (k=1; k<=nfsh;k++)
    yrs_sel_ch_fsh(k) = yrs_sel_ch_tmp(k)(1,n_sel_ch_fsh(k));
  writeinput<<"Yrs fsh_sel change: "<<yrs_sel_ch_fsh<<endl;
  for (k=1; k<=nind;k++)
    yrs_sel_ch_ind(k) = yrs_sel_ch_tmp_ind(k)(1,n_sel_ch_ind(k));
  writeinput<<"Yrs ind_sel change: "<<yrs_sel_ch_ind<<endl;
    log_sigmarprior = log(sigmarprior);
    log_input(steepnessprior);
    log_input(sigmarprior);
    nrecs_est = endyr_rec_est-styr_rec_est+1;
    nrecs_est = endyr_rec_est-styr_rec_est+1;
    writeinput<<"#  SSB estimated in styr endyr: " <<styr_sp    <<" "<<endyr_sp      <<" "<<endl;
    writeinput<<"#  Rec estimated in styr endyr: " <<styr_rec    <<" "<<endyr        <<" "<<endl;
    writeinput<<"#  SR Curve fit  in styr endyr: " <<styr_rec_est<<" "<<endyr_rec_est<<" "<<endl;
    writeinput<<"#             Model styr endyr: " <<styr        <<" "<<endyr        <<" "<<endl;
    log_qprior = log(qprior);
    writeinput<<"# qprior: " <<qprior<<" "<<endl;
    log_q_power_prior = log(q_power_prior);
    writeinput<<"# q_power_prior " <<endl<<q_power_prior<<" "<<endl;
    writeinput<<"# cv_catchbiomass " <<endl<<cv_catchbiomass<<" "<<endl;
    writeinput<<"# CatchbiomassPen " <<endl<<catchbiomass_pen<<" "<<endl;
    writeinput<<"# Number of projection years " <<endl<<nproj_yrs<<" "<<endl;// cin>>junk;

 END_CALCS
  number R_guess;

  vector offset_ind(1,nind)
  vector offset_fsh(1,nfsh)
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

  double sumtmp;
  for (k=1;k<=nfsh;k++)
    for (i=1;i<=nyrs_fsh_age(k);i++)
    {
      oac_fsh(k,i) /= sum(oac_fsh(k,i)); // Normalize to sum to one
      offset_fsh(k) -= n_sample_fsh_age(k,i)*(oac_fsh(k,i) + 0.001) * log(oac_fsh(k,i) + 0.001 ) ;
    }

  for (k=1;k<=nind;k++)
    for (i=1;i<=nyrs_ind_age(k);i++)
    {
      oac_ind(k,i) /= sum(oac_ind(k,i)); // Normalize to sum to one
      offset_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(oac_ind(k,i) + 0.001 ) ;
    }
  writeinput << "Offset index: "<<endl;
  writeinput<<offset_ind<<endl;
  writeinput << "Offset fisheries: "<<endl;
  writeinput<<offset_fsh<<endl;

  if (ad_comm::argc > 1) // Command line argument to profile Fishing mortality rates...
  {
    int on=0;
    if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-uFmort"))>-1)
      do_fmort=1;
  }

  // Compute an initial Rzero value based on exploitation 
   double btmp=0.;
   double ctmp=0.;
   dvector ntmp(1,nages);
   ntmp(1) = 1.;
   for (int a=2;a<=nages;a++)
     ntmp(a) = ntmp(a-1)*exp(-natmortprior-.05);
   btmp = wt_pop * ntmp;
   writeinput << "Mean Catch"<<endl;
   ctmp = mean(catch_bio);
   writeinput << ctmp <<endl;
   R_guess = log((ctmp/.02 )/btmp) ;
   writeinput << "R_guess "<<endl;
   writeinput << R_guess <<endl;
 END_CALCS

PARAMETER_SECTION
 // Biological Parameters
  init_bounded_number M(.02,.8,phase_M)
  matrix  natage(styr,endyr+1,1,nages)
  // vector Sp_Biom(styr_sp,endyr)
  vector pred_rec(styr_rec,endyr)
  vector mod_rec(styr_rec,endyr) // As estimated by model
  matrix  Z(styr,endyr,1,nages)
  matrix  S(styr,endyr,1,nages)
  //HM number  natmort
  vector  natmort(1,nages)
 
 // Stock rectuitment params
  init_number mean_log_rec(1); 
  init_bounded_number steepness(0.21,Steepness_UB,phase_srec)
  init_number log_Rzero(phase_Rzero)  
  // OjO
  // init_bounded_vector initage_dev(2,nages,-10,10,4)
  init_bounded_vector rec_dev(styr_rec,endyr,-10,10,2)
  // init_vector rec_dev(styr_rec,endyr,2)
  init_number log_sigmar(phase_sigmar);
  number m_sigmarsq  
  number m_sigmar
  number sigmarsq  
  number sigmar
  number alpha   
  number beta   
  number Bzero   
  number Rzero   
  number phizero
  number avg_rec_dev   

 // Fishing mortality parameters
  // init_vector         log_avg_fmort(1,nfsh,phase_fmort)
  // init_bounded_matrix fmort_dev(1,nfsh,styr,endyr,-15,15.,phase_fmort)
  init_bounded_matrix fmort(1,nfsh,styr,endyr,0.00,5.,phase_fmort)
  vector Fmort(styr,endyr);  // Annual total Fmort
  number hrate
  number catch_tmp
  number Fnew 

  !! for (k=1;k<=nfsh;k++) nselages_fsh(k)=nselages_in_fsh(k); // Sets all elements of a vector to one scalar value...
  !! for (k=1;k<=nind;k++) nselages_ind(k)=nselages_in_ind(k); // Sets all elements of a vector to one scalar value...

 //  init_3darray log_selcoffs_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_selcoff_fsh)
  !! writeinput <<"Nfsh "<< nfsh<<endl<<"Nfsh_selch "<<  n_sel_ch_fsh<<endl<<
  !! "Nselages_fsh: "<<nselages_fsh<<endl<<
  !! "PhaseSelCoffsh: "<<phase_selcoff_fsh<<endl;
  init_matrix_vector log_selcoffs_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_selcoff_fsh) // 3rd dimension out...
  init_matrix_vector  log_sel_spl_fsh(1,nfsh,1,n_sel_ch_fsh,1,nselages_fsh,phase_sel_spl_fsh)

  !! writeinput << " xx Number of fishery selectivity changes: "<<n_sel_ch_fsh<<endl;
  !! writeinput << " xx Number of index  selectivity changes: "<<n_sel_ch_ind<<endl;
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

  matrix  Ftot(styr,endyr,1,nages)
  3darray F(1,nfsh,styr,endyr,1,nages)
  3darray eac_fsh(1,nfsh,1,nyrs_fsh_age,1,nages)
  matrix  pred_catch(1,nfsh,styr,endyr)
  3darray catage(1,nfsh,styr,endyr,1,nages)
  matrix catage_tot(styr,endyr,1,nages)
  matrix expl_biom(1,nfsh,styr,endyr)

 // Parameters for computing SPR rates 
  vector F50(1,nfsh)
  vector F40(1,nfsh)
  vector F35(1,nfsh)

 // Stuff for SPR and yield projections
  number sigmar_fut
  vector ftmp(1,nfsh)
  number SB0
  number SBF50
  number SBF40
  number SBF35
  vector Fratio(1,nfsh)
  !! Fratio = 1;
  !! Fratio /= sum(Fratio);

  matrix Nspr(1,4,1,nages)
 
  matrix nage_future(styr_fut,endyr_fut,1,nages)

  init_vector rec_dev_future(styr_fut,endyr_fut,phase_proj);


 
  vector Sp_Biom_future(styr_fut-rec_age,endyr_fut);

  3darray F_future(1,nfsh,styr_fut,endyr_fut,1,nages);
  matrix Z_future(styr_fut,endyr_fut,1,nages);
  matrix S_future(styr_fut,endyr_fut,1,nages);
  matrix catage_future(styr_fut,endyr_fut,1,nages);
  number avg_rec_dev_future
  vector avg_F_future(1,5)

 // Index Observation parameters
  init_number_vector log_q_ind(1,nind,phase_q) 
  init_number_vector log_q_power_ind(1,nind,phase_q_power) 
  !! writeinput <<nind<<endl<<n_sel_ch_ind<<endl<<nselages_ind<<endl<<phase_selcoff_ind<<endl;
  init_matrix_vector log_selcoffs_ind(1,nind,1,n_sel_ch_ind,1,nselages_ind,phase_selcoff_ind)
  init_vector_vector logsel_slope_ind(1,nind,1,n_sel_ch_ind,phase_logist_ind) // Need to make positive or reparameterize
  init_vector_vector logsel_dslope_ind(1,nind,1,n_sel_ch_ind,phase_dlogist_ind) // Need to make positive or reparameterize
  matrix                sel_slope_ind(1,nind,1,n_sel_ch_ind)
  matrix                sel_dslope_ind(1,nind,1,n_sel_ch_ind)
  init_vector_vector     sel50_ind(1,nind,1,n_sel_ch_ind,phase_logist_ind)
  init_bounded_vector_vector     seld50_ind(1,nfsh,1,n_sel_ch_ind,lb_d50,nages,phase_dlogist_ind)

  3darray log_sel_ind(1,nind,styr,endyr,1,nages)
  3darray sel_ind(1,nind,styr,endyr,1,nages)
  matrix avgsel_ind(1,nind,1,n_sel_ch_ind);

  matrix pred_ind(1,nind,styr,endyr)
  3darray eac_ind(1,nind,1,nyrs_ind_age,1,nages)

 // Likelihood value names         
  number sigma
  vector rec_like(1,4)
  vector catch_like(1,nfsh)
  vector age_like_fsh(1,nfsh)
  vector age_like_ind(1,nind)
  matrix sel_like_fsh(1,nfsh,1,4)       
  matrix sel_like_ind(1,nind,1,4)       
  vector index_like(1,nind)
  vector fpen(1,6)    
  vector post_priors(1,4)
  vector post_priors_indq(1,nind)
  objective_function_value obj_fun
  vector obj_comps(1,12)

  sdreport_number B100
  sdreport_number F50_est
  sdreport_number F40_est
  sdreport_number F35_est
  vector q_ind(1,nind)
  vector q_power_ind(1,nind)
  // sdreport_vector q_ind(1,nind)
  sdreport_vector totbiom(styr,endyr+1)
  sdreport_vector totbiom_NoFish(styr,endyr)
  sdreport_vector Sp_Biom(styr_sp,endyr+1)
  sdreport_vector Sp_Biom_NoFishR(styr,endyr)
  sdreport_vector SPR(styr,endyr)
  sdreport_number ABCBiom;
  sdreport_vector recruits(styr,endyr+1)
  // vector recruits(styr,endyr+1)
  sdreport_number depletion
  sdreport_number depletion_dyn
  sdreport_number MSY;
  sdreport_number MSYL;
  sdreport_number Fmsy;
  sdreport_number lnFmsy;
  sdreport_number Fcur_Fmsy;
  sdreport_number Rmsy;
  sdreport_number Bmsy;
  sdreport_number Bcur_Bmsy;
  // NOTE TO DAVE: Need to have a phase switch for sdreport variables(
  matrix catch_future(1,4,styr_fut,endyr_fut); // Note, don't project for F=0 (it will bomb)
  sdreport_matrix SSB_fut(1,5,styr_fut,endyr_fut)
  !! writeinput <<"logRzero "<<log_Rzero<<endl;
  !! writeinput <<"logmeanrec "<<mean_log_rec<<endl;
  !! writeinput<< "exp(log_sigmarprior "<<exp(log_sigmarprior)<<endl;
 // Initialize coefficients (if needed)
 LOCAL_CALCS
  for (k=1;k<=nfsh;k++) 
  {
    writeinput<<"Fish sel phase: "<<phase_selcoff_fsh(k)<<" "<<fshname(k)<<endl;
    switch (fsh_sel_opt(k))
    {
      case 1 : // Selectivity coefficients 
      {
        if(phase_selcoff_fsh(k)>0)
        {
          writeinput<<"Initial fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            log_selcoffs_fsh(k,jj)(1,nselages_in_fsh(k)) = log_selcoffs_fsh_in(k,jj)(1,nselages_in_fsh(k));
            writeinput <<"Init coef:"<<endl<<exp(log_selcoffs_fsh(k,jj)(1,nselages_in_fsh(k))) <<endl;
          }
        }
      }
        break;
      case 2 : // Single logistic
      {
        if(phase_logist_fsh(k)<0)
        {
          writeinput<<"Fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
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
          writeinput<<"Fixing fishery sel to"<<endl<<n_sel_ch_fsh(k)<<endl;
          for (int jj=1;jj<=n_sel_ch_fsh(k);jj++) 
          {
            logsel_slope_fsh(k,jj) = logsel_slp_in_fsh(k,jj)  ;
            sel50_fsh(k,jj)        =    sel_inf_in_fsh(k,jj)  ;
          }
        }
      }
     break;
    }
  }
  for (k=1;k<=nind;k++) 
  {
    writeinput<<"Index sel phase: "<<phase_selcoff_ind(k)<<endl;
    if(phase_selcoff_ind(k)>0)
    {
      writeinput<<"Fixing "<<indname(k)<<" index sel to"<<endl<<n_sel_ch_ind(k)<<endl;
    // exit(1);
      for (int jj=1;jj<=n_sel_ch_ind(k);jj++) 
      {
        log_selcoffs_ind(k,jj)(1,nselages_in_ind(k)) = log_selcoffs_ind_in(k,jj)(1,nselages_in_ind(k));
        // writeinput <<"Init coef:"<<endl<<exp(log_selcoffs_ind(k,jj)(1,nselages_in_ind(k))) <<endl;
      }
    }
  }
  writeinput <<"Leaving data secton"<<endl;
 END_CALCS

INITIALIZATION_SECTION
  logsel_dslope_ind .9
  sel50_ind 2.
  logsel_dslope_fsh .9
  sel50_fsh 3.
  M natmortprior; 
  steepness steepnessprior
  log_sigmar log_sigmarprior;
  log_Rzero R_guess;
  mean_log_rec R_guess;
  log_q_ind log_qprior; 
  log_q_power_ind log_q_power_prior; 

  logsel_slope_fsh logsel_slp_in_fshv ;
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
  Get_Selectivity();
  Get_Mortality();
  Get_Bzero();
  Get_Numbers_at_Age();
  Get_Index_Predictions();
  Get_Fishery_Predictions();

  // Output calcs-------------------------
  if (sd_phase())
  {
    compute_spr_rates();
    if (mcmcmode)
    {
      Calc_Dependent_Vars();
      mcflag   = 0;
      mcmcmode = 0;
    }
    else
    {
      if (mcflag)
        Calc_Dependent_Vars();
    }
  }
  // Objective function calcs------------
  evaluate_the_objective_function();
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
  if (do_fmort)
    Profile_F();
 //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==

FUNCTION write_mceval
  if (mcmcmode != 3)
    write_mceval_hdr();
  mcmcmode = 3;
  get_msy();
  Future_projections();
  Calc_Dependent_Vars();
  mceval<< model_name         << " "  ;
  mceval<< obj_fun            << " "  ;
  // mceval<< sel_fsh(1,2007,12) << " "  ;
  mceval<<
  q_ind     << " "<< 
  M         << " "<< 
  steepness << " "<< 
  depletion << " "<< 
  MSY       << " "<< 
  MSYL      << " "<< 
  Fmsy      << " "<< 
  Fcur_Fmsy << " "<< 
  Bcur_Bmsy << " "<< 
  Bmsy      << " "<< 
  ABCBiom   << " "<< 
  SSB_fut(1,endyr+1)        << " "<<
  B100                      << " "<<
  SSB_fut(1,endyr+1)/B100   << " "<<
  SSB_fut(1,endyr+2)/B100   << " "<<
  SSB_fut(1,endyr+3)/B100   << " "<<
  SSB_fut(1,endyr+4)/B100   << " "<<
  SSB_fut(1,endyr+5)/B100   << " "<<
  SSB_fut(1,endyr_fut)/B100 << " "<<
  catch_future(1,endyr+1)<< " "<<
  catch_future(1,endyr+2)<< " "<<
  catch_future(1,endyr+3)<< " "<<
  catch_future(1,endyr+4)<< " "<<
  catch_future(1,endyr+5)<< " "<<
  F35       << " "<<
  F40       << " "<<
  F50       << " "<<
  SSB_fut(1,endyr_fut) << " "<< 
  SSB_fut(2,endyr_fut) << " "<< 
  SSB_fut(3,endyr_fut) << " "<< 
  SSB_fut(4,endyr_fut) << " "<< 
  SSB_fut(5,endyr_fut) << " "<< 
  catch_future(1,styr_fut)    << " "<<  
  catch_future(2,styr_fut)    << " "<<  
  catch_future(3,styr_fut)    << " "<<  
  catch_future(4,styr_fut)    << " "<<  endl;
 /*
  FUNCTION write_mceval
  if (mcmcmode != 3)
    write_mceval_hdr();
  mcmcmode = 3;
  mceval<< model_name         << " "  ;
  mceval<< obj_fun            << " "  ;
  mceval<< sel_fsh(1,2007,12) << " "  ;
  mceval<<endl;

  get_msy();
  Future_projections();
  Calc_Dependent_Vars();
  mceval<<
  q_ind     << " "<< 
  M         << " "<< 
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
  catch_future(1,styr_fut)    << " "<<  
  catch_future(2,styr_fut)    << " "<<  
  catch_future(3,styr_fut)    << " "<<  
  catch_future(4,styr_fut)    << " "<<  endl;

 */
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
        dvar_vector sel_coffs_tmp(1,nselages_fsh(k,isel_ch_tmp));
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
          log_sel_fsh(k,i)(1,nselages_fsh(k,isel_ch_tmp))        = sel_coffs_tmp;
          log_sel_fsh(k,i)(nselages_fsh(k,isel_ch_tmp),nages)    = log_sel_fsh(k,i,nselages_fsh(k,isel_ch_tmp));
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
          log_sel_fsh(k,i)(1,nselages_fsh(k,isel_ch_tmp))     = -1.*log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                                                ( age_vector(1,nselages_fsh(k,isel_ch_tmp)) - sel50_tmp) ));
          log_sel_fsh(k,i)(nselages_fsh(k,isel_ch_tmp),nages) = log_sel_fsh(k,i,nselages_fsh(k,isel_ch_tmp));
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
        log_sel_fsh(k,i)(1,nselages_fsh(k,isel_ch_tmp))     =
                     -log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                     ( age_vector(1,nselages_fsh(k,isel_ch_tmp)) - sel50_tmp) ))+
                     log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                     ( age_vector(1,nselages_fsh(k,isel_ch_tmp)) - seld50_tmp))) );

        log_sel_fsh(k,i)(nselages_fsh(k,isel_ch_tmp),nages) = 
                     log_sel_fsh(k,i,nselages_fsh(k,isel_ch_tmp));

        log_sel_fsh(k,i) -= max(log_sel_fsh(k,i));  
      }
    }
    break;
    //---Calculate the fishery selectivity from the sel_spl from nodes...
    case 4 : // Splines
    {
      /*
      int isel_ch_tmp = 1 ;
      dvar_matrix tempsel(1,n_sel_ch_fsh(k),1,nages);
      // This needs to be dimensioned by the number of changes and the number of coefficients
      writeinput <<"Splined..."<<endl; 
      vcubic_spline_function_array a(1,n_sel_ch_fsh(k),xnodes_fsh(k),log_sel_spl_fsh(k));

     tempsel =a(xages_fsh(k));

     log_sel_fsh(k,styr) = tempsel(1);;
     j=1;
     for (int i=styr+1;i<=endyr;i++) {
       if(sel_change_in_fsh(k,i)>0) {
         j++; 
         log_sel_fsh(k,i) = tempsel(j); 
       } 
       log_sel_fsh(k,i)   -= log(mean(mfexp(log_sel_fsh(k,i) )));
     }
     */
    }
    break;
    } // End of switch for fishery selectivity type
  } // End of fishery loop
  // Index specific---
  for (k=1;k<=nind;k++)
  {
    switch (ind_sel_opt(k))
    {
      case 1 : // Selectivity coefficients
      //---Calculate the fishery selectivity from the sel_coffs (Only if being used...)   
      {
        int isel_ch_tmp = 1 ;
  
        dvar_vector sel_coffs_tmp(1,nselages_ind(k,isel_ch_tmp));
        for (i=styr;i<=endyr;i++)
        {
          if (i==yrs_sel_ch_ind(k,isel_ch_tmp)) 
          {
            sel_coffs_tmp.initialize();
            sel_coffs_tmp = log_selcoffs_ind(k,isel_ch_tmp);
            avgsel_ind(k,isel_ch_tmp)              = log(mean(mfexp(sel_coffs_tmp(q_age_min(k),q_age_max(k)))));
            if (isel_ch_tmp<n_sel_ch_ind(k))
              isel_ch_tmp++;
          }
          log_sel_ind(k,i)(1,nselages_ind(k,isel_ch_tmp))        = sel_coffs_tmp;
          log_sel_ind(k,i)(nselages_ind(k,isel_ch_tmp),nages)    = log_sel_ind(k,i,nselages_ind(k,isel_ch_tmp));
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
            log_sel_ind(k,i)                                  -= log(mean(mfexp(log_sel_ind(k,i)(q_age_min(k),q_age_max(k))))); 
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
            log_sel_ind(k,i)(1,nselages_ind(k,isel_ch_tmp))     =
                         -log( 1.0 + mfexp(-1.*sel_slope_tmp * 
                         ( age_vector(1,nselages_ind(k,isel_ch_tmp)) - sel50_tmp) ))+
                         log( 1. - 1/(1.0 + mfexp(-sel_dslope_tmp * 
                         ( age_vector(1,nselages_ind(k,isel_ch_tmp)) - seld50_tmp))) );

            log_sel_ind(k,i)(nselages_ind(k,isel_ch_tmp),nages) = 
                         log_sel_ind(k,i,nselages_ind(k,isel_ch_tmp));

            log_sel_ind(k,i) -= max(log_sel_ind(k,i));  
            log_sel_ind(k,i)                                  -= log(mean(mfexp(log_sel_ind(k,i)(q_age_min(k),q_age_max(k))))); 
          }
        }
      break;
    }// end of swtiches for index selectivity
  } // End of index loop

  // Map selectivities across fisheries and indices as needed.
  for (k=1;k<=nfsh;k++)
    if (sel_map(2,k)!=k)  // If 2nd row shows a different fishery then use that fishery
      log_sel_fsh(k) = log_sel_fsh(sel_map(2,k));

  for (k=1+nfsh;k<=nfsh_and_ind;k++)
    if (sel_map(1,k)!=2) 
      log_sel_ind(k-nfsh) = log_sel_fsh(sel_map(2,k));
    else if (sel_map(2,k)!=(k-nfsh)) 
      log_sel_ind(k-nfsh) = log_sel_ind(sel_map(2,k));

  sel_fsh = mfexp(log_sel_fsh);
  sel_ind = mfexp(log_sel_ind);

FUNCTION Get_Mortality
  natmort = M + M_offset;
  if (!Popes)
  {
    Fmort.initialize();
    Z.initialize();
    for (i=styr;i<=endyr;i++)
    {
      for (k=1;k<=nfsh;k++)
      {
        Fmort(i) +=  fmort(k,i);
        F(k,i)    =  fmort(k,i) * sel_fsh(k,i) ;
        Z(i)    += F(k,i);
      }
      Z(i) += natmort;
    }
    S  = mfexp(-1.*Z);
  }
  

FUNCTION Get_Numbers_at_Age
  // natage(styr,1) = mfexp(mean_log_rec + rec_dev(styr)); 
  // Recruitment in subsequent years
  for (i=styr+1;i<=endyr;i++)
    natage(i,1)=mfexp(mean_log_rec+rec_dev(i));

  mod_rec(styr)  = natage(styr,1);

  if (Popes)
  {
    dvar_vector  t1=mfexp(-natmort*0.5);
    dvar_vector  t2=mfexp(-natmort);
    for (i=styr;i<=endyr;i++)
    {
      Catch_at_Age(i);
      // Pope's approximation //   Next year N     =   This year x NatSurvivl - catch
      natage(i+1)(2,nages) = ++(elem_prod(natage(i)(1,nages-1),t2(1,nages-1)) - elem_prod(catage_tot(i)(1,nages-1),t1(1,nages-1)));
      Ftot(i)(1,nages-1) = log(natage(i)(1,nages-1)) - --log(natage(i+1)(2,nages)) - natmort(2,nages);
      natage(i+1,nages)   += natage(i,nages)*t2(nages) - catage_tot(i,nages)*t1(nages);
      // Approximation to "F" continuous form for computing within-year sp biomass
      Ftot(i,nages)      = log(natage(i,nages-1)+natage(i,nages)) -log(natage(i+1,nages)) -natmort(nages);
      // writeinput <<i<<" "<<Ftot(i)(nages-4,nages)<<endl; // cout <<i<<" "<<natage(i)<<endl; // cout <<i<<" "<<natage(i+1)<<endl;
      dvariable ctmp=sum(catage_tot(i));
      for (k=1;k<=nfsh;k++)
      {
        F(k,i)  = Ftot(i) * sum(catage(k,i))/ctmp;
      }
      Z(i)    = Ftot(i)+natmort;
      S(i)    = mfexp(-Z(i));
    }
  }
  else // Baranov
  {
    for (i=styr;i<=endyr;i++)
    {
      // get_Fs( i ); //ojo, add switch here for different catch equation XX
      if (i!=endyr)
      {
        natage(i+1)(2,nages) = ++elem_prod(natage(i)(1,nages-1),S(i)(1,nages-1));
        natage(i+1,nages)+=natage(i,nages)*S(i,nages);
      }
      Catch_at_Age(i);
      Sp_Biom(i)  = elem_prod(natage(i),pow(S(i),spmo_frac)) * wt_mature; 
      // if (current_phase()>2)cout <<(i)<<" "<<spmo_frac<<" "<<wt_mature<<" "<<S(i)<<endl;
      if (i<endyr) mod_rec(i+1)  = natage(i+1,1);
    }
  }

FUNCTION Get_Index_Predictions
  // Index computations------------------
  dvariable sum_tmp;
  sum_tmp.initialize();
  int yy;
  for (k=1;k<=nind;k++)
  {
    for (i=styr;i<=endyr;i++)
    {        
      pred_ind(k,i) = q_ind(k) * pow(elem_prod(natage(i),pow(S(i),ind_month_frac(k))) * 
                                     elem_prod(sel_ind(k,i) , wt_ind(k,i)),q_power_ind(k));
    }
    for (i=1;i<=nyrs_ind_age(k);i++)
    {
      yy = int(yrs_ind_age(k,i)); 
      dvar_vector tmp_n   = elem_prod(pow(S(yy),ind_month_frac(k)),elem_prod(sel_ind(k,yy),natage(yy)));  
      sum_tmp             = sum(tmp_n);
      if (use_age_err)
        eac_ind(k,i)      = age_err * tmp_n/sum_tmp;
      else
        eac_ind(k,i)      = tmp_n/sum_tmp;

    }
  }
FUNCTION Get_Fishery_Predictions
  for (k=1;k<=nfsh;k++)
  {
    for (i=1; i<=nyrs_fsh_age(k); i++)
      if (use_age_err)
        eac_fsh(k,i) = age_err * catage(k,yrs_fsh_age(k,i))/sum(catage(k,yrs_fsh_age(k,i)));
      else
        eac_fsh(k,i) = catage(k,yrs_fsh_age(k,i))/sum(catage(k,yrs_fsh_age(k,i)));
  }

FUNCTION Calc_Dependent_Vars
  // cout<<"In DepVar stage 1"<<endl;
  get_msy();
  if (phase_proj>0) 
    Future_projections();

  dvar_matrix N_nofsh(styr,endyr,1,nages);
  N_nofsh.initialize();
  N_nofsh(styr) = natage(styr);
  for (i=styr;i<=endyr;i++)
  {                 
    recruits(i)  = natage(i,1);
    if (i>styr)
    {
      N_nofsh(i,1) = recruits(i);
      N_nofsh(i)(2,nages) = ++(elem_prod(N_nofsh(i-1)(1,nages-1),exp(-natmort(1,nages-1))));
      N_nofsh(i,nages)   += N_nofsh(i-1,nages)*exp(-natmort(nages));
    }
    totbiom_NoFish(i) = N_nofsh(i)*wt_pop;
    totbiom(i)        = natage(i)*wt_pop;
    Sp_Biom_NoFishR(i) = Sp_Biom(i) / (elem_prod(N_nofsh(i),pow(exp(-natmort),spmo_frac)) * wt_mature); 
    depletion         = totbiom(endyr)/totbiom(styr);
    depletion_dyn     = totbiom(endyr)/totbiom_NoFish(endyr);
    dvar_vector F_agetmp(1,nages);
    F_agetmp.initialize();
    for (k=1;k<=nfsh;k++)
      F_agetmp += F(k,i);
    SPR(i) = Implied_SPR(F_agetmp);
  }
  B100 = phizero * mean(recruits(styr_rec_est, endyr_rec_est));
  dvar_vector Ntmp(1,nages);
  Ntmp(2,nages) = ++elem_prod(natage(endyr)(1,nages-1),S(endyr)(1,nages-1));
  Ntmp(nages)  += natage(endyr,nages)*S(endyr,nages);
  Ntmp(1)       = SRecruit(Sp_Biom(endyr+1-rec_age));
  ABCBiom       = Ntmp*wt_pop;
  recruits(endyr+1) = Ntmp(1);
  totbiom(endyr+1)  = ABCBiom;

FUNCTION void Catch_at_Age(const int& i)
  dvariable vbio=0.;
  dvariable pentmp;
  dvar_vector Nmid(1,nages);
  dvar_vector Ctmp(1,nages);
  catage_tot(i).initialize();
  if (Popes)
  {
    // Nmid = natage(i)*mfexp(-natmort/2. ); 
    Nmid = natage(i)*mfexp(-M/2. ); 
  }
  for (k=1;k<=nfsh;k++)
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
      catage_tot(i) += Ctmp;                      
      catage(k,i)    = Ctmp;                      
      if (last_phase())
        pred_catch(k,i) = Ctmp*wt_fsh(k,i);
    }
    else
    {
      catage(k,i) = elem_prod(elem_div(F(k,i),Z(i)),elem_prod(1.-S(i),natage(i)));
      pred_catch(k,i) = catage(k,i)*wt_fsh(k,i);
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
  if (active(log_Rzero)) // OjO
    obj_fun += .5 * square(log_Rzero-mean_log_rec); // A slight penalty to keep Rzero in reality...

  obj_comps.initialize();
  obj_comps(1) = sum(catch_like);
  obj_comps(2) = sum(age_like_fsh);
  obj_comps(3) = sum(sel_like_fsh);
  obj_comps(4) = sum(index_like);
  obj_comps(5) = sum(age_like_ind);
  obj_comps(6) = sum(sel_like_ind);
  obj_comps(7) = sum(rec_like);
  obj_comps(8) = sum(fpen);
  obj_comps(9) = sum(post_priors_indq);
  obj_comps(10)= sum(post_priors);
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
  if (active(rec_dev))
  {
    sigmar     =  mfexp(log_sigmar);
    sigmarsq   =  square(sigmar);
    dvariable SSQRec;
    SSQRec.initialize();
    if (current_phase()>2)
    {
      if (last_phase())
      {
        pred_rec = SRecruit(Sp_Biom(styr_rec-rec_age,endyr-rec_age).shift(styr_rec)(styr_rec,endyr));
      }
      else 
      {
         pred_rec = 1.+SRecruit(Sp_Biom(styr_rec-rec_age,endyr-rec_age).shift(styr_rec)(styr_rec,endyr));
      }
      dvar_vector chi(styr_rec_est,endyr_rec_est);
      chi = log(mod_rec(styr_rec_est,endyr_rec_est)) - log(pred_rec(styr_rec_est,endyr_rec_est));
      // cout<<(mod_rec(styr_rec_est,endyr_rec_est)) <<endl<< (pred_rec(styr_rec_est,endyr_rec_est)+.001)<<endl << (Sp_Biom(styr_rec-rec_age,endyr-rec_age).shift(styr_rec)(styr_rec,endyr))<<endl; ;exit(1);
      SSQRec   = norm2( chi ) ;
      m_sigmar   =  sqrt( SSQRec  / nrecs_est);
      m_sigmarsq =  m_sigmar * m_sigmar   ;

      if (current_phase()>4||last_phase())
        rec_like(1) = (SSQRec+ sigmarsq/2.)/(2*sigmarsq) + nrecs_est*log_sigmar; 
      else
        rec_like(1) = .1*(SSQRec+ sigmarsq/2.)/(2*sigmarsq) + nrecs_est*log_sigmar; 
    }

    if (last_phase())
    {
      // Variance term for the parts not estimated by sr curve
      rec_like(4) += .5*norm2( rec_dev(styr_rec,styr_rec_est) )/sigmarsq + (styr_rec_est-styr_rec)*log(sigmar) ; 

      if ( endyr > endyr_rec_est)
        rec_like(4) += .5*norm2( rec_dev(endyr_rec_est,endyr  ) )/sigmarsq + (endyr-endyr_rec_est)*log(sigmar) ; 
    }
    else // JNI comment next line
      rec_like(2) += norm2( rec_dev(styr_rec_est,endyr) ) ;

    // rec_like(2) += norm2( rec_dev(styr_rec_est,endyr) ) ;

    if (active(rec_dev_future))
    {
      // Future recruitment variability (based on past)
      sigmar_fut   = sigmar ;
      rec_like(3) += norm2(rec_dev_future)/(2*square(sigmar_fut))+ size_count(rec_dev_future)*log(sigmar_fut);
    }
    // cout <<rec_like<<" "<< SSQRec<<" "<<sigmarsq<<endl;
  }

FUNCTION Compute_priors
  post_priors.initialize();
  post_priors_indq.initialize();
  for (k=1;k<=nind;k++)
  {
    if (active(log_q_ind(k)))
      post_priors_indq(k) += square(q_ind(k)-qprior(k))/(2*cvqprior(k)*cvqprior(k)); 
    if (active(log_q_power_ind(k)))
      post_priors_indq(k) += square(q_power_ind(k)-q_power_prior(k))/(2*cvq_power_prior(k)*cvq_power_prior(k)); 
  }

  if (active(M))
    post_priors(1) += square(M-natmortprior)/(2*cvnatmortprior*cvnatmortprior); 

  if (active(steepness))
    post_priors(2) += square(steepness-steepnessprior)/(2*cvsteepnessprior*cvsteepnessprior); 

  if (active(log_sigmar))
    post_priors(3) += square(sigmar-sigmarprior)/(2*cvsigmarprior*cvsigmarprior); 

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
    if (active(log_selcoffs_fsh(k)))
    {
      for (i=1;i<=n_sel_ch_fsh(k);i++)
      {
        int iyr = yrs_sel_ch_fsh(k,i) ;
        sel_like_fsh(k,1) += curv_pen_fsh(k)*norm2(first_difference(
                                                   first_difference(log_sel_fsh(k,iyr ))));
        // This part is the penalty on the change itself--------------
        if (i>1)
        {
          dvariable var_tmp = square(sel_change_in_fsh(k,iyr ));
          sel_like_fsh(k,2)    += .5*norm2( log_sel_fsh(k,iyr-1) - log_sel_fsh(k,iyr) ) / var_tmp ;
        }
        int nagestmp = nselages_fsh(k,1);
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
    if (active(log_selcoffs_ind(k)))
    {
      for (i=1;i<=n_sel_ch_ind(k);i++)
      {
        int iyr = yrs_sel_ch_ind(k,i) ;
        sel_like_ind(k,1) += curv_pen_ind(k)*norm2(first_difference(
                                                 first_difference(log_sel_ind(k,iyr))));
        // This part is the penalty on the change itself--------------
        if (i>1)
        {
          dvariable var_tmp = square(sel_change_in_ind(k,iyr ));
          sel_like_ind(k,2)    += .5*norm2( log_sel_ind(k,iyr-1) - log_sel_ind(k,iyr) ) 
                                   / var_tmp ;
        }
        int nagestmp = nselages_ind(k,1);
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
  index_like.initialize();
  int yy;
  for (k=1;k<=nind;k++)
    for (i=1;i<=nyrs_ind(k);i++)
    {
      yy = int(yrs_ind(k,i));
      index_like(k) += square(log(obs_ind(k,i)) - log(pred_ind(k,yy)) ) / 
                                   (2.*obs_lse_ind(k,i)*obs_lse_ind(k,i));
    }
  /* normal distribution option to add someday...
    for (i=1;i<=nyrs_ind(k);i++)
      index_like(k) += square(obs_ind(k,i) - pred_ind(k,yrs_ind(k,i)) ) / 
                                   (2.*obs_se_ind(k,i)*obs_se_ind(k,i));
  */

FUNCTION Age_Like
  age_like_fsh.initialize();
  for (k=1;k<=nfsh;k++)
    for (int i=1;i<=nyrs_fsh_age(k);i++)
      age_like_fsh(k) -= n_sample_fsh_age(k,i)*(oac_fsh(k,i) + 0.001) * log(eac_fsh(k,i) + 0.001 ) ;
  age_like_fsh -= offset_fsh;

  age_like_ind.initialize();
  for (k=1;k<=nind;k++)
    for (int i=1;i<=nyrs_ind_age(k);i++)
      age_like_ind(k) -= n_sample_ind_age(k,i)*(oac_ind(k,i) + 0.001) * log(eac_ind(k,i) + 0.001 ) ;
  age_like_ind -= offset_ind;

FUNCTION Oper_Model
 // Initialize things used here only
  Calc_Dependent_Vars();
  Write_Datafile();

  int nsims;
  ifstream sim_in("nsims.dat");
  sim_in >> nsims; sim_in.close();

  dvector ran_ind_vect(1,nind);
  ofstream SaveOM("Om_Out.dat",ios::app);
  double C_tmp;
  dvariable Fnow;
  dvariable meanrec;
  meanrec=mean(recruits);
  dvector new_ind(1,nind);
  new_ind.initialize();
  system("cls"); cout<<"Number of replicates: "<<endl;
  // Initialize recruitment in first year
  nage_future(styr_fut,1) = meanrec;
  nage_future(styr_fut)(2,nages)              = ++elem_prod(natage(endyr)(1,nages-1),S(endyr)(1,nages-1));
  nage_future(styr_fut,nages)                += natage(endyr,nages)*S(endyr,nages);
  for (int isim=1;isim<=nsims;isim++)
  {
    cout<<isim<<" ";
    // Copy file to get mean for Mgt Strategies
    system("init_stuff.bat");
    for (i=styr_fut;i<=endyr_fut;i++)
    {
      // Some unit normals...
      ran_ind_vect.fill_randn(rng);
      // Create new index observations
      for (k = 1 ; k<= nind ; k++)
        new_ind(k) = mfexp(ran_ind_vect(k)*.2)*value(nage_future(i)*q_ind(k)*sel_ind(k,endyr)); // use value function since converts to a double
      // Append new index observation to datafile
      ofstream ind_out("NewSrv.dat",ios::app);
      ind_out <<i<<" "<< new_ind<<endl; ind_out.close();
      system("ComputeTAC.bat"); // commandline function to get TAC (catchnext.dat)

      // Now read in TAC (actual catch)
      ifstream CatchNext("CatchNext.dat");
      CatchNext >> C_tmp; CatchNext.close();
      Fnow = SolveF2(endyr,nage_future(i), C_tmp);

      F_future(1,i) = sel_fsh(1,endyr) * Fnow;
      Z_future(i)   = F_future(1,i) + natmort;
      S_future(i)   = mfexp(-Z_future(i));
      // nage_future(i,1)  = SRecruit( Sp_Biom_future(i-rec_age) ) * mfexp(rec_dev_future(i)) ;     

      nage_future(i,1)  = meanrec;
      Sp_Biom_future(i) = wt_mature * elem_prod(nage_future(i),pow(S_future(i),spmo_frac)) ;
      // Now graduate for the next year....
      if (i<endyr_fut)
      {
        nage_future(i+1)(2,nages) = ++elem_prod(nage_future(i)(1,nages-1),S_future(i)(1,nages-1));
        nage_future(i+1,nages)   += nage_future(i,nages)*S_future(i,nages);
      }
      catage_future(i) = 0.; 
      for (k = 1 ; k<= nfsh ; k++)
        catage_future(i) += elem_prod(nage_future(i) , elem_prod(F_future(k,i) , elem_div( ( 1.- S_future(i) ) , Z_future(i))));
  
      SaveOM << model_name       <<
        " "  << isim             <<
        " "  << i                <<
        " "  << Fnow             <<
        " "  << nage_future(i)(5,nages)*wt_fsh(1,endyr)(5,nages)<<
        " "  << Sp_Biom_future(i)<<
        " "  << catage_future(i)*wt_fsh(1,endyr)<<
        " "  << nage_future(i)<<
        " "  << F_future(1,i)<<
      endl;
    }
  }
  SaveOM.close();
  if (!mceval_phase())
    exit(1);

FUNCTION void get_future_Fs(const int& i,const int& iscenario)
    // get F's
    switch (iscenario)
    {
      case 1:
        ftmp = F50;
        // cout <<"F35 "<<i<<" "<<ftmp<<endl;
        break;
      case 2:
        ftmp = F40;
        // ftmp = SolveF2(endyr,nage_future(i), 0.1  * sum(catch_bioT(endyr)) );
        break;
      case 3:
        ftmp = F35;
        // ftmp = SolveF2(endyr,nage_future(i), 0.05 * sum(catch_bioT(endyr)) );
        break;
      case 4:
        // ftmp = SolveF2(endyr,nage_future(i), 0.01 * sum(catch_bioT(endyr)) );
        ftmp = Fmsy;
        break;
      case 5:
        ftmp = 0.0;
        break;
    }
    Z_future(i) = natmort;
    for (k=1;k<=nfsh;k++)
    {
      F_future(k,i) = sel_fsh(k,endyr) * ftmp(k);
      Z_future(i)  += F_future(k,i);
    }
    S_future(i) = mfexp(-Z_future(i));

FUNCTION Future_projections
  // Need to check on treatment of Fratio--whether it should be included or not
  SSB_fut.initialize();
  catch_future.initialize();
  for (int iscen=1;iscen<=5;iscen++)
  {
   // Future Sp_Biom set equal to estimated Sp_Biom w/ right lag
    // Sp_Biom_future(styr_fut-rec_age,styr_fut-1) = Sp_Biom(endyr-rec_age+1,endyr);
    for (i=styr_fut-rec_age;i<=styr_fut-1;i++)
      Sp_Biom_future(i) = wt_mature * elem_prod(natage(i),pow(S(i),spmo_frac)) ;

    // cout<<Sp_Biom(endyr-10,endyr)<<endl<<Sp_Biom_future<<endl;exit(1);
    nage_future(styr_fut)(2,nages) = ++elem_prod(natage(endyr)(1,nages-1),S(endyr)(1,nages-1));
    nage_future(styr_fut,nages)   += natage(endyr,nages)*S(endyr,nages);
    Sp_Biom_future(styr_fut)       = wt_mature * elem_prod(nage_future(i),pow(S_future(i),spmo_frac)) ;
    // Future Recruitment (and Sp_Biom)
    for (i=styr_fut;i<endyr_fut;i++)
    {
      nage_future(i,1)  = SRecruit( Sp_Biom_future(i-rec_age) ) * mfexp(rec_dev_future(i)) ;     
      get_future_Fs(i,iscen);
      // Now graduate for the next year....
      nage_future(i+1)(2,nages) = ++elem_prod(nage_future(i)(1,nages-1),S_future(i)(1,nages-1));
      nage_future(i+1,nages)   += nage_future(i,nages)*S_future(i,nages);
      Sp_Biom_future(i) = wt_mature * elem_prod(nage_future(i),pow(S_future(i),spmo_frac)) ;
    }
    nage_future(endyr_fut,1)  = SRecruit( Sp_Biom_future(endyr_fut-rec_age) ) * mfexp(rec_dev_future(endyr_fut)) ;     
    get_future_Fs(endyr_fut,iscen);
    Sp_Biom_future(endyr_fut)  = wt_mature * elem_prod(nage_future(endyr_fut),pow(S_future(endyr_fut),spmo_frac)) ;
    // Now get catch at future ages
    for (i=styr_fut; i<=endyr_fut; i++)
    {
      catage_future(i) = 0.; 
      for (k = 1 ; k<= nfsh ; k++)
      {
        catage_future(i) += elem_prod(nage_future(i) , elem_prod(F_future(k,i) , 
                              elem_div( ( 1.- S_future(i) ) , Z_future(i))));
        if (iscen!=5) 
          catch_future(iscen,i)   += catage_future(i)*wt_fsh(k,endyr);
      }
      SSB_fut(iscen,i) = Sp_Biom_future(i);
      // cout<<iscen<<" "<<i<<" "<< SSB_fut(iscen,i) <<" "<<nage_future(i)(1,4)<<endl;
    }
  }   //End of loop over F's
  Sp_Biom(endyr+1) = Sp_Biom_future(styr_fut);

FUNCTION get_msy
 /*Function calculates used in calculating MSY and MSYL for a designated component of the
  population, given values for stock recruitment and selectivity...  
  Fmsy is the trial value of MSY example of the use of "funnel" to reduce the amount of storage for derivative calculations */

  dvariable sumF=0.;
  for (k=1;k<=nfsh;k++)
    sumF += sum(F(k,endyr));
  for (k=1;k<=nfsh;k++)
    Fratio(k) = sum(F(k,endyr)) / sumF;

  dvariable Stmp;
  dvariable Rtmp;
  double df=1.e-05;
  dvariable F1;
  F1.initialize();
  F1 = (1.1*natmortprior);
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
    yld1   = yield(Fratio,F1);
    yld2   = yield(Fratio,F2);
    yld3   = yield(Fratio,F3);
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
    ttt      = yld(Fratio,F1);
    Fmsy     = F1;
    MSY      = ttt(2);
    Bmsy     = ttt(1);
    MSYL     = ttt(1)/Bzero;
    lnFmsy   = log(MSY/ttt(5)); // Exploitation fraction relative to total biomass
    Bcur_Bmsy= Sp_Biom(endyr)/Bmsy;

    dvariable FFtmp;
    FFtmp.initialize();
    for (k=1;k<=nfsh;k++)
      FFtmp += mean(F(k,endyr));
    Fcur_Fmsy= FFtmp/Fmsy;
    Rmsy     = Rtmp;
  }

 //+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+ 
FUNCTION dvar_vector yld(const dvar_vector& Fratio, const dvariable& Ftmp)
  RETURN_ARRAYS_INCREMENT();
  /*dvariable utmp=1.-mfexp(-(Ftmp)); dvariable Ntmp; dvariable Btmp; dvariable yield; dvariable survtmp=exp(-1.*natmort); dvar_vector seltmp=sel_fsh(endyr); Ntmp = 1.; Btmp = Ntmp*wt(1)*seltmp(1); Stmp = .5*Ntmp*wt(1)*maturity(1); yield= 0.; for ( j=1 ; j < nages ; j++ ) { Ntmp  *= (1.-utmp*seltmp(j))*survtmp; Btmp  += Ntmp*wt(j+1)*seltmp(j+1); Stmp  += .5 * Ntmp *wt(j+1)*maturity(j+1); } //Max Age - 1 yr yield   += utmp * Btmp; Ntmp    /= (1-survtmp*(1.-utmp*seltmp(nages))); Btmp    += Ntmp*wt(nages)*seltmp(nages); Stmp    += 0.5 *wt(nages)* Ntmp *maturity(nages); yield   += utmp * Btmp; //cout<<yield<<" "<<Stmp<<" "<<Btmp<<" ";*/
  dvar_vector msy_stuff(1,5);
  dvariable phi;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  msy_stuff.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  for (k=1;k<=nfsh;k++)
   seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = natmort;
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
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) ) * wt_mature;
  // Req    = Requil(phi) * exp(sigmarsq/2);
  msy_stuff(5)  = Ntmp * wt_pop;      
  msy_stuff(4)  = phi/phizero ;       // SPR
  msy_stuff(3)  = Requil(phi) ;       // Eq Recruitment
  msy_stuff(5) *= msy_stuff(3);       // BmsyTot
  msy_stuff(2) *= msy_stuff(3);       // MSY
  msy_stuff(1)  = phi*(msy_stuff(3)); // Bmsy
  RETURN_ARRAYS_DECREMENT();
  return msy_stuff;

FUNCTION dvariable yield(const dvar_vector& Fratio, const dvariable& Ftmp)
  RETURN_ARRAYS_INCREMENT();
  /*dvariable utmp=1.-mfexp(-(Ftmp)); dvariable Ntmp; dvariable Btmp; dvariable yield; dvariable survtmp=exp(-1.*natmort); dvar_vector seltmp=sel_fsh(endyr); Ntmp = 1.; Btmp = Ntmp*wt(1)*seltmp(1); Stmp = .5*Ntmp*wt(1)*maturity(1); yield= 0.; for ( j=1 ; j < nages ; j++ ) { Ntmp  *= (1.-utmp*seltmp(j))*survtmp; Btmp  += Ntmp*wt(j+1)*seltmp(j+1); Stmp  += .5 * Ntmp *wt(j+1)*maturity(j+1); } //Max Age - 1 yr yield   += utmp * Btmp; Ntmp    /= (1-survtmp*(1.-utmp*seltmp(nages))); Btmp    += Ntmp*wt(nages)*seltmp(nages); Stmp    += 0.5 *wt(nages)* Ntmp *maturity(nages); yield   += utmp * Btmp; //cout<<yield<<" "<<Stmp<<" "<<Btmp<<" ";*/
  dvariable phi;
  dvariable Req;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  dvariable   yield;
  yield.initialize();

  dvar_matrix seltmp(1,nfsh,1,nages);
  for (k=1;k<=nfsh;k++)
   seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = natmort;
  for (k=1;k<=nfsh;k++)
  { 
    Fatmp(k) = Fratio(k) * Ftmp * seltmp(k);
    Ztmp    += Fatmp(k);
  } 
  dvar_vector survtmp = mfexp(-Ztmp);
  // cout<<Ftmp<<" ";

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
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) )* wt_mature;
  // Req    = Requil(phi) * mfexp(sigmarsq/2);
  Req    = Requil(phi) ;
  yield *= Req;

  RETURN_ARRAYS_DECREMENT();
  return yield;

FUNCTION dvariable yield(const dvar_vector& Fratio, dvariable& Ftmp, dvariable& Stmp,dvariable& Req)
  RETURN_ARRAYS_INCREMENT();
  dvariable phi;
  dvar_vector Ntmp(1,nages);
  dvar_vector Ctmp(1,nages);
  dvariable   yield   = 0.;

  dvar_matrix seltmp(1,nfsh,1,nages);
  for (k=1;k<=nfsh;k++)
   seltmp(k) = sel_fsh(k,endyr); // NOTE uses last-year of fishery selectivity for projections.

  dvar_matrix Fatmp(1,nfsh,1,nages);
  dvar_vector Ztmp(1,nages);

  Ztmp = natmort;
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
  phi    = elem_prod( Ntmp , pow(survtmp,spmo_frac ) )* wt_mature;
  // Req    = Requil(phi) * exp(sigmarsq/2);
  Req    = Requil(phi) ;
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
 dvariable sumF=0.;
  for (k=1;k<=nfsh;k++)
    sumF += sum(F(k,endyr));
  for (k=1;k<=nfsh;k++)
    Fratio(k) = sum(F(k,endyr)) / sumF;
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
  prof_F <<0.0<<" "<< Bzero <<" "<<0.0<<" "<<Rzero<< " 1.00"<<endl; 
  dvar_vector ttt(1,5);
  for (int ii=1;ii<=500;ii++)
  {
    F1    = double(ii)/500;
    yld1  = yield(Fratio,F1,Stmp,Rtmp);
    ttt   = yld(Fratio,F1);
    prof_F <<F1<<" "<< ttt << endl; 
  } 

FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
  switch (SrType)
  {
    case 1:
      RecTmp = elem_prod((Stmp / phizero) , mfexp( alpha * ( 1. - Stmp / Bzero ))) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = elem_prod(Stmp , 1. / ( alpha + beta * Stmp));        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec);                    //Avg recruitment
      break;
    case 4:
      RecTmp = elem_prod(Stmp , mfexp( alpha  - Stmp * beta)) ; //Old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(const double& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero) * mfexp( alpha * ( 1. - Stmp / Bzero )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha + beta * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec);                    //Avg recruitment
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha  - Stmp * beta) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(const dvariable& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero) * mfexp( alpha * ( 1. - Stmp / Bzero )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha + beta * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(mean_log_rec );                    //Avg recruitment
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha  - Stmp * beta) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

 //=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
FUNCTION Get_Bzero
  Bzero.initialize();
  Rzero    =  mfexp(log_Rzero); 

  dvar_vector survtmp = mfexp(-natmort);

  dvar_matrix natagetmp(styr_rec,styr,1,nages);
  natagetmp.initialize();

  natagetmp(styr_rec,1) = Rzero;
  for (j=2; j<=nages; j++)
  {
    natagetmp(styr_rec,j) = natagetmp(styr_rec,j-1) * survtmp(j-1);
    Bzero += wt_mature(j) * pow(survtmp(j),spmo_frac)*natagetmp(styr_rec,j) ;
  }
  natagetmp(styr_rec,nages) /= (1.-survtmp(nages)); 
  Bzero += wt_mature(nages) * pow(survtmp(nages),spmo_frac)*natagetmp(styr_rec,nages) ;

  phizero = Bzero/Rzero;

  switch (SrType)
  {
    case 1:
      alpha = log(-4.*steepness/(steepness-1.));
      break;
    case 2:
    {
      alpha  =  Bzero * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero;
      beta   = (5. * steepness - 1.) / (4. * steepness * Rzero);
    }
    break;
    case 4:
    {
      beta  = log(5.*steepness)/(0.8*Bzero) ;
      alpha = log(Rzero/Bzero)+beta*Bzero;
    }
      break;
  }
  Sp_Biom.initialize();
  Sp_Biom(styr_sp,styr_rec) = Bzero;
  for (i=styr_rec;i<styr;i++)
  {
    natagetmp(i,1) = mfexp(rec_dev(i) + mean_log_rec);
    for (j=1; j<=nages; j++)
      Sp_Biom(i) += natagetmp(i,j)*pow(survtmp(j),spmo_frac) * wt_mature(j); 

    natagetmp(i+1)(2,nages) = ++(elem_prod(natagetmp(i)(1,nages-1),mfexp(-natmort(1,nages-1))));
    natagetmp(i+1,nages)   += natagetmp(i,nages)*mfexp(-natmort(nages));
  }
  natagetmp(styr,1)   = mfexp(rec_dev(styr) + mean_log_rec);
  mod_rec(styr_rec,styr) = column(natagetmp,1);
  natage(styr)  = natagetmp(styr); // OjO
  for (j=1; j<=nages; j++)
    Sp_Biom(styr) += natagetmp(styr,j)*pow(survtmp(j),spmo_frac) * wt_mature(j); 
  // Sp_Biom(styr) = natagetmp(styr)*pow(surv,spmo_frac) * wt_mature; 
  // cout <<natagetmp<<endl;exit(1);

FUNCTION dvariable Requil(dvariable& phi)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp =  Bzero * (alpha + log(phi) - log(phizero) ) / (alpha*phi);
      break;
    case 2:
      RecTmp =  (phi-alpha)/(beta*phi);
      break;
    case 3:
      RecTmp =  mfexp(mean_log_rec);
      break;
    case 4:
      RecTmp =  (log(phi)+alpha) / (beta*phi); //RecTmp =  (log(phi)/alpha + 1.)*beta/phi;
      break;
  }
  // Req    = Requil(phi) * exp(sigmarsq/2);
  // return RecTmp* exp(sigmarsq/2);
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION write_mceval_hdr
    mceval<< "Model Postererior ";
    for (k=1;k<=nind;k++)
      mceval<< " q_ind_"<< k<< " ";
    mceval<<"M steepness depletion MSY MSYL Fmsy Fcur_Fmsy Bcur_Bmsy Bmsy totbiom_"<<endyr<<" "<< 
    " SSB_yr_"<<endyr+1 <<" "<<
    " B100         "<< 
    " B_"<<endyr+1   <<"_over_B100 " << 
    " B_"<<endyr+2   <<"_over_B100 " << 
    " B_"<<endyr+3   <<"_over_B100 " << 
    " B_"<<endyr+4   <<"_over_B100 " << 
    " B_"<<endyr+5   <<"_over_B100 " << 
    " B_"<<endyr_fut <<"_over_B100 " << 
    " fut_catch_"<<endyr+1<<" "<<  
    " fut_catch_"<<endyr+2<<" "<<  
    " fut_catch_"<<endyr+3<<" "<<  
    " fut_catch_"<<endyr+4<<" "<<  
    " fut_catch_"<<endyr+5<<" "<<  
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
    for (k=1;k<=nind;k++)
    {
      cout<<indname(k)<<endl;
      cout<<get_AC(k)<<endl<<endl;
    }
    if (!Popes)
      for (k=1;k<=nfsh;k++)
        Ftot += F(k);
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
        Fmort(i) += mean(F(k,i));
    report << Fmort<<endl;
    report << "Total mortality (Z)"<<endl;
    report << Z<<endl;
    report << "Estimated numbers of fish " << endl;
    for (i=styr;i<=endyr;i++) 
      report <<"       Year: "<< i << " "<< natage(i) << endl;
    report << endl<< "Estimated F mortality " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "Fishery "<< k <<" : "<< endl ;
      for (i=styr;i<=endyr;i++) 
        report << "        Year: "<<i<<" "<<F(k,i)<<  " "<< endl;
    }

    report << endl<< "Observed index values " << endl;
    for (k=1;k<=nind;k++)
    {
      int ii=1;
      report <<endl<< "Yr_Obs_Pred_Index "<< k <<" : "<< endl ;
      for (i=styr;i<=endyr;i++)
      {
        if (ii<=yrs_ind(k).indexmax())
        {
          if (yrs_ind(k,ii)==i)
          {
            report << i<< " "<< obs_ind(k,ii) << " "<< pred_ind(k,i) <<endl;
            ii++;
          }
          else
            report << i<< " -1 "<< " "<< pred_ind(k,i)<<endl;
        }
        else
          report << i<< " -1 "<< " "<< pred_ind(k,i)<<endl;
      }
    }

    report << endl<< "Index_Q:  "<<q_ind << endl;

    report << endl<< "Observed Prop " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "ObsFishery "<< k <<" : "<< endl ;
      for (i=1;i<=nyrs_fsh_age(k);i++) 
        report << yrs_fsh_age(k,i)<< " "<< oac_fsh(k,i) << endl;
    }
    report << endl<< "Predicted prop  " << endl;
    for (k=1;k<=nfsh;k++)
    {
      report << "PredFishery "<< k <<" : "<< endl;
      for (i=1;i<=nyrs_fsh_age(k);i++) 
        report << yrs_fsh_age(k,i)<< " "<< eac_fsh(k,i) << endl;
    }
    report << endl<< "Observed prop Index" << endl;
    for (k=1;k<=nind;k++)
    {
      report << "ObsIndex "<<k<<" : "<<  endl;
      for (i=1;i<=nyrs_ind_age(k);i++) 
        report << yrs_ind_age(k,i)<< " "<< oac_ind(k,i) << endl;
    }
    report << endl<< "Predicted prop Index" << endl;
    for (k=1;k<=nind;k++)
    {
      report << "PredIndex "<<k<<" : "<<  endl;
      for (i=1;i<=nyrs_ind_age(k);i++) 
        report << yrs_ind_age(k,i)<< " "<< eac_ind(k,i) << endl;
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
        report << "Index  "<< k <<"  "<< i<<" "<<sel_ind(k,i) << endl;

    report << endl<< "Stock Recruitment stuff "<< endl;
    for (i=styr_rec;i<=endyr;i++)
      if (active(log_Rzero))
        report << i<< " "<<Sp_Biom(i-rec_age)<< " "<< SRecruit(Sp_Biom(i-rec_age))<< " "<< mod_rec(i)<<endl;
      else 
        report << i<< " "<<Sp_Biom(i-rec_age)<< " "<< " 999" << " "<< mod_rec(i)<<endl;

    report << endl<< "Curve to plot "<< endl;
    report <<"stock Recruitment"<<endl;
    report <<"0 0 "<<endl;
    dvariable stock;
    for (i=1;i<=30;i++)
    {
      stock = double (i) * Bzero /25.;
      if (active(log_Rzero))
        report << stock <<" "<< SRecruit(stock)<<endl;
      else
        report << stock <<" 99 "<<endl;
    }

    report   << endl<<"Likelihood Components" <<endl;
    report   << "----------------------------------------- " <<endl;
    report   << "  catch_like  age_like_fsh sel_like_fsh index_like age_like_ind sel_like_ind rec_like fpen post_priors_indq post_priors residual total"<<endl;
    report   << " "<<obj_comps<<endl;

    obj_comps(11)= obj_fun - sum(obj_comps) ; // Residual 
    obj_comps(12)= obj_fun ;                  // Total
    report   <<"  catch_like       "<<setw(10)<<obj_comps(1) <<endl
             <<"  age_like_fsh     "<<setw(10)<<obj_comps(2) <<endl
             <<"  sel_like_fsh     "<<setw(10)<<obj_comps(3) <<endl
             <<"  index_like        "<<setw(10)<<obj_comps(4) <<endl
             <<"  age_like_ind     "<<setw(10)<<obj_comps(5) <<endl
             <<"  sel_like_ind     "<<setw(10)<<obj_comps(6) <<endl
             <<"  rec_like         "<<setw(10)<<obj_comps(7) <<endl
             <<"  fpen             "<<setw(10)<<obj_comps(8) <<endl
             <<"  post_priors_indq "<<setw(10)<<obj_comps(9) <<endl
             <<"  post_priors      "<<setw(10)<<obj_comps(10)<<endl
             <<"  residual         "<<setw(10)<<obj_comps(11)<<endl
             <<"  total            "<<setw(10)<<obj_comps(12)<<endl;

    report   << endl;
    report << "Fit to Catch Biomass "<<endl;
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
  
    report   << "index Likelihood(s) " <<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nind;k++)
      report << "  Index_#"<< k <<"  " << index_like(k)<<endl;
    report   << endl;

    report << setw(10)<< setfixed() << setprecision(5) <<endl;
    report   << "Age likelihoods for indices :"<<endl;
    report   << "-------------------------" <<endl;
    for (k=1;k<=nind;k++)
      report << "  Age_Index_#"<< k <<"  " << age_like_ind(k)<<endl;
    report   << endl;

    report   << "Selectivity penalties for indices :"<<endl;
    report   << "-------------------------" <<endl;
    report   << "  Index Curvature_Age Change_Time Dome_Shaped"<<endl;
    for (k=1;k<=nind;k++)
      report << "  Sel_Index_#"<< k <<"  "<< sel_like_ind(k,1) <<" "<<sel_like_ind(k,2)<<" "<<sel_like_ind(k,3)<< endl;
    report   << endl;

    report << setw(10)<< setfixed() << setprecision(5) <<endl;
    report   << "Recruitment penalties: " <<rec_like<<endl;
    report   << "-------------------------" <<endl;
    report   << "  (sigmar)            " <<sigmar<<endl;
    report   << "  S-R_Curve           " <<rec_like(1)<< endl;
    report   << "  Regularity          " <<rec_like(2)<< endl;
    report   << "  Future_Recruits     " <<rec_like(3)<< endl;
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
      report << "Q_Index_#"<< k <<"           "
             << setw(10)<<post_priors_indq(k) 
             << setw(10)<< q_ind(k)
             << setw(10)<< qprior(k)
             << setw(10)<< cvqprior(k)<<endl;

      report << "Q_power_Index_#"<< k <<"           "
             << setw(10)<<post_priors_indq(k) 
             << setw(10)<< q_power_ind(k)
             << setw(10)<< q_power_prior(k)
             << setw(10)<< cvq_power_prior(k)<<endl;
    }

    // writerep(post_priors(1),repstring);
    // cout <<repstring<<endl;
    report   << "Natural_Mortality     "
             << setw(10)<< post_priors(1)
             << setw(10)<< M
             << setw(10)<< natmortprior
             << setw(10)<< cvnatmortprior <<endl;
    report   << "Steepness             "
             << setw(10)<< post_priors(2)
             << setw(10)<< steepness
             << setw(10)<< steepnessprior
             << setw(10)<< cvsteepnessprior <<endl;
    report   << "SigmaR                "
             << setw(10)<< post_priors(3)
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
  report<<"SR_Curve_fit__in_styr_endyr: " <<styr_rec_est<<" "<<endyr_rec_est<<" "<<endl;
  report<<"Model_styr_endyr:            " <<styr        <<" "<<endyr        <<" "<<endl;

  report<<"M_prior,_CV,_phase "<< natmortprior<< " "<< cvnatmortprior<<" "<<phase_M<<endl;
  report<<"qprior,_CV,_phase " <<qprior<<" "<<cvqprior<<" "<< phase_q<<endl;
  report<<"q_power_prior,_CV,_phase " <<q_power_prior<<" "<<cvq_power_prior<<" "<< phase_q_power<<endl;

  report<<"cv_catchbiomass: " <<cv_catchbiomass<<" "<<endl;
  report<<"Projection_years "<< nproj_yrs<<endl;
  for (k=1;k<=nfsh;k++)
    report << "Fsh_sel_opt_fish: "<<k<<" "<<fsh_sel_opt(k)<<" "<<sel_change_in_fsh(k)<<endl;
  for (k=1;k<=nind;k++)
    report<<"Index_Sel_Opt_: " <<k<<" "<<(ind_sel_opt(k))<<endl;
    
  report <<"Phase_index_Sel_Coffs: "<<phase_selcoff_ind<<endl; 
  report <<"Fshry_Selages: " << nselages_in_fsh  <<endl;
  report <<"Index_Selages: " << nselages_in_ind <<endl;
  report << "Phase_for_age-spec_fishery "<<phase_selcoff_fsh<<endl;
  report << "Phase_for_logistic_fishery "<<phase_logist_fsh<<endl;
  report << "Phase_for_dble_logistic_fishery "<<phase_dlogist_fsh<<endl;
  report << "Phase_for_age-spec_index  "<<phase_selcoff_ind<<endl;
  report << "Phase_for_logistic_index  "<<phase_logist_ind<<endl;
  report << "Phase_for_dble_logistic_indy "<<phase_dlogist_ind<<endl;

  for (k=1; k<=nfsh;k++)
  {
    report <<"Number_of_select_changes_fishery: "<<k<<" "<<n_sel_ch_fsh(k)<<endl;
    report<<"Yrs_fsh_sel_change: "<<yrs_sel_ch_fsh(k)<<endl;
    report << "sel_change_in: "<<sel_change_in_fsh(k) << endl;
  }
  for (k=1; k<=nind;k++)
  {
    report <<"Number_of_select_changes_index: "<<k<<" "<<n_sel_ch_ind(k)<<endl;
    report<<"Yrs_ind_sel_change: "<<yrs_sel_ch_ind(k)<<endl;
    report << "sel_change_in: "<<sel_change_in_ind(k) << endl;
  }

FUNCTION write_msy_out
  ofstream msyout("msyout.dat");
  msyout << " # Natural Mortality       " <<endl;
  for (j=1;j<=nages;j++) 
    msyout <<M <<" ";
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
  projout << nfsh<<"   # Number of fisheries"<<endl;
  projout <<"14   # Number of projection years"<<endl;
  projout <<"1000 # Number of simulations"<<endl;
  projout <<endyr<< " # Begin year of projection" <<endl;
  projout <<nages<< " # Number of ages" <<endl;
  for (j=1;j<=nages;j++) 
    projout <<M <<" ";
  projout << " # Natural Mortality       " <<endl;
  double sumtmp;
  sumtmp = 0.;
  for (k=1;k<=nfsh;k++) 
    sumtmp += catch_bio(k,endyr);
  projout << sumtmp<< " # TAC in current year (assumed catch) " <<endl;
  projout << sumtmp<< " # TAC in current year+1 (assumed catch) " <<endl;

  for (k=1;k<=nfsh;k++) 
    projout <<  F(k,endyr)/mean((F(k,endyr)))<<" ";
   //  + fmort_dev(k,endyr)) /Fmort(endyr)<<" ";

  projout << "   # Fratio                  " <<endl;

  projout << mean(Fmort(endyr-4,endyr))<<"  # average f                  " <<endl;
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
    projout<< sel_fsh(k,endyr) <<" ";
  projout<< endl;
  projout <<"# natage"<<endl<< natage(endyr) << endl;
  // projout <<"#_N_recruitment_years (not including last 3 estimates)"<<endl<<endyr-(1977+rec_age+3) << endl;
  // projout <<"#_Recruitment_start_at_1977_yearclass=1978_for_age_1_recruits"<<yy(1977+rec_age,endyr-3)<<endl<<mod_rec(1977+rec_age,endyr-3)<< endl;

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
 dvector seltmp = value(sel_fsh(1,endyr)) 
                 +value(sel_fsh(1,endyr-1))  
                 +value(sel_fsh(1,endyr-2))  
                 +value(sel_fsh(1,endyr-3))  
                 +value(sel_fsh(1,endyr-4));  
 newproj << mean(Fmort(endyr-4,endyr))<<endl;
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
   for (j=1;j<=nages;j++) newproj <<natmort(j) <<" "; newproj<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<maturity<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt_pop<< endl;
 newproj <<"#_Wt_at_age_fishery" <<aa<<endl<<wt_fsh(1,endyr) << endl;
 newproj <<"#" <<endl;

 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl;
 seltmp = value(sel_fsh(1,endyr)) +value(sel_fsh(1,endyr-1))  +value(sel_fsh(1,endyr-2));  
 newproj << seltmp/max(seltmp)<<endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage(endyr)<< endl;
 // newproj <<"#_N_recruitment_years (not including last estimate)"<<endl<<endyr-(1977+rec_age) << endl;

 // newproj <<"#_Recruitment_start_at_1977_yearclass=1978_for_age_1_recruits"<<yy(1977+rec_age,endyr-1) <<endl<<mod_rec(1977+rec_age,endyr-1)<< endl;

 newproj <<"#_Spawning biomass "<<endl<<Sp_Biom(styr-rec_age,endyr-rec_age)/1000<< endl;
 newproj.close();
 
RUNTIME_SECTION
  convergence_criteria 1.e-1,,1.e-01,1.e-03,1e-5,1e-5
  maximum_function_evaluations 1500

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1500);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1500);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize=500000000;

FINAL_SECTION
  // Calc_Dependent_Vars();
  Write_R();
  write_proj();
  write_projout();
  write_msy_out();
  Profile_F();
  
GLOBALS_SECTION
  #include <admodel.h>
	/**
	\def log_input(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef log_input
	#define log_input(object) writeinput << "# " #object "\n" << object << endl;
	
  ofstream writeinput("writeinput.log");
 // void get_sel_changes(int& k);
  adstring_array fshname;
  adstring_array indname;
  adstring simname;
  adstring model_name;
  adstring projfile_name;
  adstring datafile_name;
  adstring cntrlfile_name;
  adstring tmpstring;
  adstring repstring;
  adstring version_info;

FUNCTION dvariable Implied_SPR( const dvar_vector& F_age) 
  // Function that returns SPR percentage given a realized value of F...
    dvar_vector ntmp0(1,nages);
    dvar_vector ntmp1(1,nages);
    dvar_vector survivtmp=mfexp(-natmort);
    ntmp0(1) = 1.;
    ntmp1(1) = 1.;
    for (j=2;j<nages;j++)
    {
      ntmp0(j)  = ntmp0(j-1)* survivtmp(j-1);
      ntmp1(j)  = ntmp1(j-1)* exp(-1.*(natmort(j-1) + F_age(j-1) ));
    }
    ntmp0(nages)  =  ntmp0(nages-1)* survivtmp(nages)/ (1.- survivtmp(nages));
    ntmp1(nages)  =  ntmp1(nages-1)* exp(-(natmort(nages) + F_age(nages-1)))/ (1.- mfexp(-(natmort(nages) + F_age(nages) )));
    dvariable sb0_tmp;
    dvariable sb1_tmp;

    sb0_tmp.initialize();
    sb1_tmp.initialize();
    for (j=1;j<=nages;j++)
    {
      // natmort till spawning 
      sb0_tmp += ntmp0(j)*wt_mature(j) * mfexp(-spmo_frac * natmort(j));
      sb1_tmp += ntmp1(j)*wt_mature(j) * mfexp(-spmo_frac * ( natmort(j) + F_age(j) ));
    }
    return(sb1_tmp / sb0_tmp);
 
FUNCTION dvariable get_spr_rates(double spr_percent)
  RETURN_ARRAYS_INCREMENT();
  dvar_matrix sel_tmp(1,nages,1,nfsh);
  sel_tmp.initialize();
  for (k=1;k<=nfsh;k++)
    for (j=1;j<=nages;j++)
      sel_tmp(j,k) = sel_fsh(k,endyr,j); // NOTE uses last-year of fishery selectivity for projections.
  dvariable sumF=0.;
  for (k=1;k<=nfsh;k++)
  {
    Fratio(k) = sum(F(k,endyr)) ;
    sumF += Fratio(k) ;
  }
  Fratio /= sumF;
  double df=1.e-3;
  dvariable F1 ;
  F1.initialize();
  F1 = .8*natmortprior;
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
    yld1   = -1000*square(log(spr_percent/spr_ratio(F1, sel_tmp)));
    yld2   = -1000*square(log(spr_percent/spr_ratio(F2, sel_tmp)));
    yld3   = -1000*square(log(spr_percent/spr_ratio(F3, sel_tmp)));
    dyld   = (yld2 - yld3)/(2*df);                          // First derivative (to find the root of this)
    dyldp  = (yld3-(2*yld1)+yld2)/(df*df);  // Newton-Raphson approximation for second derivitive
    F1    -= dyld/dyldp;
  }
  RETURN_ARRAYS_DECREMENT();
  return(F1);

FUNCTION dvariable spr_ratio(dvariable trial_F,dvar_matrix sel_tmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable SBtmp;
  dvar_vector Ntmp(1,nages);
  dvar_vector indtmp(1,nages);
  SBtmp.initialize();
  Ntmp.initialize();
  indtmp.initialize();
  dvar_matrix Ftmp(1,nages,1,nfsh);
  Ftmp = sel_tmp;
  for (j=1;j<=nages;j++) 
  {
    Ftmp(j) = elem_prod(Ftmp(j), trial_F * Fratio);
    indtmp(j)  = mfexp(-sum(Ftmp(j)) - natmort(j));
  }
  Ntmp(1)=1.;
  j=1;
  SBtmp  += Ntmp(j)*wt_mature(j)*pow(indtmp(j),spmo_frac);
  for (j=2;j<nages;j++)
  {
    Ntmp(j) = Ntmp(j-1)*indtmp(j-1);
    SBtmp  += Ntmp(j)*wt_mature(j)*pow(indtmp(j),spmo_frac);
  }
  Ntmp(nages)=Ntmp(nages-1)*indtmp(nages-1)/(1.-indtmp(nages));

  SBtmp  += Ntmp(nages)*wt_mature(nages)*pow(indtmp(nages),spmo_frac);
  RETURN_ARRAYS_DECREMENT();
  return(SBtmp/phizero);

FUNCTION dvariable spr_unfished()
  RETURN_ARRAYS_INCREMENT();
  dvariable Ntmp;
  dvariable SBtmp;
  SBtmp.initialize();
  Ntmp = 1.;
  for (j=1;j<nages;j++)
  {
    SBtmp += Ntmp*wt_mature(j)*exp(-spmo_frac * natmort(j));
    Ntmp  *= mfexp( -natmort(j));
  }
  Ntmp    /= (1.-exp(-natmort(nages)));
  SBtmp += Ntmp*wt_mature(j)*exp(-spmo_frac * natmort(nages));
  RETURN_ARRAYS_DECREMENT();
  return(SBtmp);

FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  dvariable sumF=0.;
  for (k=1;k<=nfsh;k++)
  {
    Fratio(k) = sum(F(k,endyr)) ;
    sumF += Fratio(k) ;
  }
  Fratio /= sumF;

  F35_est = get_spr_rates(.35);
  F50_est = get_spr_rates(.50);
  F40_est = get_spr_rates(.40);

  for (k=1;k<=nfsh;k++)
  {
    F50(k) = F50_est * (Fratio(k));
    F40(k) = F40_est * (Fratio(k));
    F35(k) = F35_est * (Fratio(k));
  }

  /* FUNCTION get_agematrix
  // B3+($B$5-$B$4)/(COUNT($B$2:$L$2)-1)
  for (i=1;i<=nages;i++)
    {
     sd_age(i)=cv_age(i)*len_age(i);
     var_age(i)=sd_age(i)*sd_age(i);
     for (int j=1;j<=nlenbins;j++)
       {
        diff = len(j)-len_age(i);
        trans(i,j)=2/sd_age(i)*exp(-diff*diff/(2*var_age(i)));
       }
      trans(i)=trans(i)/sum(trans(i));
    }
    */ 


FUNCTION void writerep(dvariable& tmp,adstring& tmpstring)
  cout <<tmpstring<<endl<<endl;
  tmpstring = printf("3.5%f",value(tmp));

FUNCTION void do_Newton_Raphson_for_mortality(dvariable hrate)
  dvariable Fold ;
  Fold = hrate;
  for (int ii=1;ii<=4;ii++)
  {
      dvariable ZZ = Fold + M;
      dvariable SS = mfexp(-ZZ);
      dvariable AA = Fold * (1. - SS);
      dvariable BB = ZZ;
      dvariable CC = 1. + (Fold - 1) * SS;
      dvariable dd = 1.;
      dvariable FX = AA / BB - hrate;
      dvariable FPX = (BB * CC - AA * dd) / (BB * BB);
      Fnew = Fold - FX / FPX;
      Fold = Fnew;
  }

FUNCTION dvariable SolveF2(const int& iyr, const dvar_vector& N_tmp, const double&  TACin)
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc; 
  dvar_matrix Fratsel(1,nfsh,1,nages);
  dvar_vector Z_tmp(1,nages) ;
  dvar_vector S_tmp(1,nages) ;
  dvar_vector Ftottmp(1,nages);
  dvariable btmp =  N_tmp * elem_prod(sel_fsh(1,iyr),wt_pop);
  dvariable ftmp;
  ftmp = TACin/btmp;
    for (k=1;k<=nfsh;k++)
      Fratsel(k) = Fratio(k)*sel_fsh(k,iyr);
    for (int ii=1;ii<=5;ii++)
    {
      Ftottmp.initialize();
      for (k=1;k<=nfsh;k++)
        Ftottmp += ftmp*Fratsel(k);
  
      Z_tmp = Ftottmp  + M; 
      S_tmp = mfexp( -Z_tmp );
      cc = 0.0;
      for (k=1;k<=nfsh;k++)
        cc += wt_fsh(k,endyr) * elem_prod(elem_div(ftmp*Fratsel(k),  Z_tmp),elem_prod(1.-S_tmp,N_tmp)); // Catch equation (vectors)
  
      dd = cc / TACin - 1.;
      if (dd<0.) dd *= -1.;
      ftmp += (TACin-cc) / btmp;
    }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp);

FUNCTION dvar_vector SolveF2(const int& iyr, const dvector&  Catch)
  // Returns vector of F's (given year) by fleet
  // Requires: N and fleet specific wts & selectivities at age, catch 
  // iterate to get Z's right
  RETURN_ARRAYS_INCREMENT();
  dvariable dd = 10.;
  dvariable cc; 
  dvar_matrix  seltmp(1,nfsh,1,nages);
  dvar_matrix  wt_tmp(1,nfsh,1,nages);
  dvar_matrix Fratsel(1,nfsh,1,nages);
  dvar_vector N_tmp = natage(iyr);
  dvar_vector Z_tmp(1,nages) ;
  dvar_vector S_tmp(1,nages) ;
  dvar_vector Ftottmp(1,nages);
  dvar_vector Frat(1,nfsh);
  dvar_vector btmp(1,nfsh);
  dvar_vector ftmp(1,nfsh);
  dvar_vector hrate(1,nfsh);
  btmp.initialize(); 
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
    for (k=1;k<=nfsh;k++)
    {
      if (hrate(k) <.9999) 
      {
        for (int ii=1;ii<=6;ii++)
        {
          Ftottmp.initialize();
          Ftottmp   = ftmp*Fratsel;
          Z_tmp = Ftottmp  + M; 
          S_tmp = mfexp( -Z_tmp );
          cc = wt_tmp(k) * elem_prod(elem_div(ftmp(k)*Fratsel(k),  Z_tmp),elem_prod(1.-S_tmp,N_tmp)); // Catch equation (vectors)
          ftmp(k) += ( Catch(k)-cc ) / btmp(k);
        }
        Frat(k)    = ftmp(k)/sum(ftmp);
        Fratsel(k) = Frat(k)*seltmp(k);
      }
    }
  }
  RETURN_ARRAYS_DECREMENT();
  return(ftmp);


FUNCTION Write_Datafile
  int nsims;
  // get the number of simulated datasets to create...
  ifstream sim_in("nsims.dat"); sim_in >> nsims; sim_in.close();
  char buffer [33];
  // compute the autocorrelation term for residuals of fit to indices...
  for (k=1;k<=nind;k++)
    ac(k) = get_AC(k);
  for (int isim=1;isim<=nsims;isim++)
  {
    // Create the name of the simulated dataset
    // simname = "sim_"+ adstring(itoa(isim,buffer,10)) + ".dat";
    simname = "sim_crap.dat";
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
        p = value(catage(k,iyr));
        p /= sum(p);
        bin.fill_multinomial(rng,p); // fill a vector v
        for (int j=1;j<=n_sample_fsh_age(k,i);j++)
          freq(bin(j))++;
        // Apply ageing error to samples..............
        p = age_err *freq/sum(freq); 
        cout << p  <<endl;
        simdat << p  <<endl;
        // Compute total catch given this sample size
        Ctmp = catch_bio(k,iyr) / (p*wt_fsh(k,iyr)); 
        // Simulated catage = proportion sampled
        sim_catage(k,i) = p * Ctmp;
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
    dmatrix new_ind(1,nind,1,nyrs_ind);
    new_ind.initialize();
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
                     value(elem_prod(wt_ind(k,iyr),elem_prod(pow(S(iyr),ind_month_frac(k)), natage(iyr)))*
                     q_ind(k)*sel_ind(k,iyr)); 
      // do next years correlated with previous
      for (i=2;i<=nyrs_ind(k);i++)
      {
        iyr=yrs_ind(k,i);
        corr_dev(k,i) = ac(k) * corr_dev(k,i-1) + sqrt(1.-square(ac(k))) * corr_dev(k,i);
        new_ind(k,i) = mfexp(corr_dev(k,i) * obs_lse_ind(k,i) ) * 
                        value(elem_prod(wt_ind(k,iyr),elem_prod(pow(S(iyr),ind_month_frac(k)), 
                        natage(iyr))) * q_ind(k)*sel_ind(k,iyr)); 
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
        p = age_err * value(elem_prod( elem_prod(pow(S(iyr),ind_month_frac(k)), natage(iyr))*q_ind(k) , sel_ind(k,iyr))); 
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
        avail_biom(i) = wt_fsh(k,i)*value(elem_prod(natage(i),sel_fsh(k,i))); 
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
        simdat<<fshname(k)<<" "<<yrs_fsh_age(k,i)<<" "<<sim_catage(k,i) <<endl;
    }
  }
  exit(1);
  // End of simulating datasets...................
  

FUNCTION Write_R
  ofstream R_report("For_R.dat");
  R_report<<"$Yr"<<endl; for (i=styr;i<=endyr;i++) R_report<<i<<" "; R_report<<endl;

  R_report<<"$TotF"<<endl << Ftot<<endl;
  R_report<<"$natmort"<<endl << natmort<<endl;

  R_report<<"$TotBiom_NoFish"<<endl; for (i=styr;i<=endyr;i++) 
  {
    double lb=value(totbiom_NoFish(i)/exp(2.*sqrt(log(1+square(totbiom_NoFish.sd(i))/square(totbiom_NoFish(i))))));
    double ub=value(totbiom_NoFish(i)*exp(2.*sqrt(log(1+square(totbiom_NoFish.sd(i))/square(totbiom_NoFish(i))))));
    R_report<<i<<" "<<totbiom_NoFish(i)<<" "<<totbiom_NoFish.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }
  R_report<<"$SSB_NoFishR"<<endl; for (i=styr;i<=endyr;i++) 
  {
    double lb=value(Sp_Biom_NoFishR(i)/exp(2.*sqrt(log(1+square(Sp_Biom_NoFishR.sd(i))/square(Sp_Biom_NoFishR(i))))));
    double ub=value(Sp_Biom_NoFishR(i)*exp(2.*sqrt(log(1+square(Sp_Biom_NoFishR.sd(i))/square(Sp_Biom_NoFishR(i))))));
    R_report<<i<<" "<<Sp_Biom_NoFishR(i)<<" "<<Sp_Biom_NoFishR.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }

  R_report<<"$TotBiom"<<endl; 
  for (i=styr;i<=endyr;i++) 
  {
    double lb=value(totbiom(i)/exp(2.*sqrt(log(1+square(totbiom.sd(i))/square(totbiom(i))))));
    double ub=value(totbiom(i)*exp(2.*sqrt(log(1+square(totbiom.sd(i))/square(totbiom(i))))));
    R_report<<i<<" "<<totbiom(i)<<" "<<totbiom.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }

  for (int k=1;k<=5;k++){
    R_report<<"$SSB_fut_"<<k<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      double lb=value(SSB_fut(k,i)/exp(2.*sqrt(log(1+square(SSB_fut.sd(k,i))/square(SSB_fut(k,i))))));
      double ub=value(SSB_fut(k,i)*exp(2.*sqrt(log(1+square(SSB_fut.sd(k,i))/square(SSB_fut(k,i))))));
      R_report<<i<<" "<<SSB_fut(k,i)<<" "<<SSB_fut.sd(k,i)<<" "<<lb<<" "<<ub<<endl;
    }
  }
  double ctmp;
  for (int k=1;k<=5;k++){
    R_report<<"$Catch_fut_"<<k<<endl; 
    for (i=styr_fut;i<=endyr_fut;i++) 
    {
      if (k==5) ctmp=0.;else ctmp=value(catch_future(k,i));
      R_report<<i<<" "<<ctmp<<endl;
    }
  }

  R_report<<"$SSB"<<endl; for (i=styr_sp;i<=endyr+1;i++) 
  {
    double lb=value(Sp_Biom(i)/exp(2.*sqrt(log(1+square(Sp_Biom.sd(i))/square(Sp_Biom(i))))));
    double ub=value(Sp_Biom(i)*exp(2.*sqrt(log(1+square(Sp_Biom.sd(i))/square(Sp_Biom(i))))));
    R_report<<i<<" "<<Sp_Biom(i)<<" "<<Sp_Biom.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }

  R_report<<"$R"<<endl; for (i=styr;i<=endyr;i++) 
  {
    double lb=value(recruits(i)/exp(2.*sqrt(log(1+square(recruits.sd(i))/square(recruits(i))))));
    double ub=value(recruits(i)*exp(2.*sqrt(log(1+square(recruits.sd(i))/square(recruits(i))))));
    R_report<<i<<" "<<recruits(i)<<" "<<recruits.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }
    R_report << "$N"<<endl;
    for (i=styr;i<=endyr;i++) 
      R_report <<   i << " "<< natage(i) << endl;
      R_report   << endl;

    for (int k=1;k<=nfsh;k++)
    {
      R_report << "$F_age_"<< (k) <<""<< endl ;
      for (i=styr;i<=endyr;i++) 
        R_report <<i<<" "<<F(k,i)<<" "<< endl;
        R_report   << endl;
    }

    R_report <<endl<< "$Fshry_names"<< endl;
    for (int k=1;k<=nfsh;k++)
      R_report << fshname(k) << endl ;

    R_report <<endl<< "$Index_names"<< endl;
    for (int k=1;k<=nind;k++)
      R_report << indname(k) << endl ;

    for (int k=1;k<=nind;k++)
    {
      int ii=1;
      R_report <<endl<< "$Obs_Index_"<< k <<""<< endl ;
      for (i=styr;i<=endyr;i++)
      {
        if (ii<=yrs_ind(k).indexmax())
        {
          if (yrs_ind(k,ii)==i)
          {
            R_report << i<< " "<< obs_ind(k,ii) << " "<< pred_ind(k,i) <<" "<< obs_se_ind(k,ii) <<endl; //values of index value (annual)
            ii++;
          }
          else
            R_report << i<< " -1 "<< " "<< pred_ind(k,i)<<" -1 "<<endl;
        }
        else
          R_report << i<< " -1 "<< " "<< pred_ind(k,i)<<" -1 "<<endl;
      }
      R_report   << endl;
    }

    R_report << endl<< "$Index_Q"<<endl;
    R_report<< q_ind << endl;
    R_report   << endl;
    for (int k=1;k<=nfsh;k++)
    {
      if (nyrs_fsh_age(k)>0) 
      { 
        R_report << "$pobs_fsh_"<< (k) <<""<< endl;
        for (i=1;i<=nyrs_fsh_age(k);i++) 
          R_report << yrs_fsh_age(k,i)<< " "<< oac_fsh(k,i) << endl;
        R_report   << endl;
      }
    }
    for (int k=1;k<=nfsh;k++)
    {
      if (nyrs_fsh_age(k)>0) 
      { 
        R_report << "$phat_fsh_"<< (k) <<""<< endl;
        for (i=1;i<=nyrs_fsh_age(k);i++) 
          R_report << yrs_fsh_age(k,i)<< " "<< eac_fsh(k,i) << endl;
          R_report   << endl;
      }
    }
    for (int k=1;k<=nind;k++)
    {
      if (nyrs_ind_age(k)>0) 
      { 
        R_report << "$pobs_ind_"<<(k)<<""<<  endl;
        for (i=1;i<=nyrs_ind_age(k);i++) 
          R_report << yrs_ind_age(k,i)<< " "<< oac_ind(k,i) << endl;
          R_report   << endl;
      }
    }
    for (int k=1;k<=nind;k++)
    {
      if (nyrs_ind_age(k)>0) 
      { 
        R_report << "$phat_ind_"<<(k)<<""<<  endl;
        for (i=1;i<=nyrs_ind_age(k);i++) 
          R_report << yrs_ind_age(k,i)<< " "<< eac_ind(k,i) << endl;
          R_report   << endl;
      }
    }
    for (int k=1;k<=nfsh;k++)
    {
      R_report << endl<< "$Obs_catch_"<<(k) << endl;
      R_report << catch_bio(k) << endl;
      R_report   << endl;
      R_report << "$Pred_catch_" <<(k) << endl;
      R_report << pred_catch(k) << endl;
      R_report   << endl;
    }

    for (int k=1;k<=nfsh;k++)
    {
      R_report << "$F_fsh_"<<(k)<<" "<<endl;
      for (i=styr;i<=endyr;i++)
      {
        R_report<< i<< " ";
        R_report<< mean(F(k,i)) <<" "<< mean(F(k,i))*max(sel_fsh(k,i)) << " ";
        R_report<< endl;
      }
    }

    for (int k=1;k<=nfsh;k++)
    {
      R_report << endl<< "$sel_fsh_"<<(k)<<"" << endl;
      for (i=styr;i<=endyr;i++)
        R_report << k <<"  "<< i<<" "<<sel_fsh(k,i) << endl; 
      R_report   << endl;
    }

    for (int k=1;k<=nind;k++)
    {
      R_report << endl<< "$sel_ind_"<<(k)<<"" << endl;
      for (i=styr;i<=endyr;i++)
        R_report << k <<"  "<< i<<" "<<sel_ind(k,i) << endl;
        R_report << endl;

    }
    R_report << endl<< "$Stock_Rec"<< endl;
    for (i=styr_rec;i<=endyr;i++)
      if (active(log_Rzero))
        R_report << i<< " "<<Sp_Biom(i-rec_age)<< " "<< SRecruit(Sp_Biom(i-rec_age))<< " "<< mod_rec(i)<<endl;
      else 
        R_report << i<< " "<<Sp_Biom(i-rec_age)<< " "<< " 999" << " "<< mod_rec(i)<<endl;
        
        R_report   << endl;

    R_report <<"$stock_Rec_Curve"<<endl;
    R_report <<"0 0"<<endl;
    dvariable stock;
    for (i=1;i<=30;i++)
    {
      stock = double (i) * Bzero /25.;
      if (active(log_Rzero))
        R_report << stock <<" "<< SRecruit(stock)<<endl;
      else
        R_report << stock <<" 99 "<<endl;
    }
    R_report   << endl;

    R_report   << endl<<"$Like_Comp" <<endl;
    obj_comps(11)= obj_fun - sum(obj_comps) ; // Residual
    obj_comps(12)= obj_fun ;

    R_report   <<obj_comps<<endl;
    R_report   << endl;
    R_report   << endl<<"$Like_Comp_names" <<endl;
    R_report   <<"catch_like     "<<endl
             <<"age_like_fsh     "<<endl
             <<"sel_like_fsh     "<<endl
             <<"index_like        "<<endl
             <<"age_like_ind     "<<endl
             <<"sel_like_ind     "<<endl
             <<"rec_like         "<<endl
             <<"fpen             "<<endl
             <<"post_priors_indq "<<endl
             <<"post_priors      "<<endl
             <<"residual         "<<endl
             <<"total            "<<endl;
    for (int k=1;k<=nfsh;k++)
    {
      R_report << "$Sel_Fshry_"<< (k) <<""<<endl;
      R_report << sel_like_fsh(k) << endl;
    }
    R_report   << endl;
  
    for (int k=1;k<=nind;k++)
    {
      R_report << "$Index_"<< (k) <<"" <<endl;
      R_report<< index_like(k)<<endl;
    }
    R_report   << endl;

    R_report << setw(10)<< setfixed() << setprecision(5) <<endl;
    for (int k=1;k<=nind;k++)
    {
      R_report << "$Age_Index_"<< (k) <<"" <<endl;
      R_report << age_like_ind(k)<<endl;
    }
    R_report   << endl;

    for (int k=1;k<=nind;k++)
    {
      R_report << "$Sel_Index_"<< (k) <<""<<endl;
      R_report<< sel_like_ind(k,1) <<" "<<sel_like_ind(k,2)<<" "<<sel_like_ind(k,3)<< endl;
    }
    R_report   << endl;

    R_report << setw(10)<< setfixed() << setprecision(5) <<endl;
    R_report   << "$Rec_Pen" <<endl<<sigmar<<"  "<<rec_like<<endl;
    R_report   << endl;

    R_report   << "$F_Pen" <<endl;
    R_report<<fpen(1)<<"  "<<fpen(2)<<endl;
    R_report   << endl;
    for (int k=1;k<=nind;k++)
    {
      R_report << "$Q_Index_"<< (k) <<""<<endl
             << " "<<post_priors_indq(k)
             << " "<< q_ind(k)
             << " "<< qprior(k)
             << " "<< cvqprior(k)<<endl;
      R_report << "$Q_power_Index_"<< (k) <<""<<endl
             << " "<<post_priors_indq(k)
             << " "<< q_power_ind(k)
             << " "<< q_power_prior(k)
             << " "<< cvq_power_prior(k)<<endl;
    }
             R_report   << endl;
    R_report << "$M_prior"<<endl;
    R_report << " "<< post_priors(1)
             << " "<< M
             << " "<< natmortprior
             << " "<< cvnatmortprior <<endl;
    R_report   << endl;
    R_report << "$Steep"<<endl;
    R_report << " "<< post_priors(2)
             << " "<< steepness
             << " "<< steepnessprior
             << " "<< cvsteepnessprior <<endl;
    R_report   << endl;
    R_report << "$Sigmar"<<endl;
    R_report << " "<< post_priors(3)
             << " "<< sigmar
             << " "<< sigmarprior
             << " "<< cvsigmarprior <<endl;
    R_report   << endl;
    R_report<<"$Num_parameters_Est"<<endl;
    R_report<<initial_params::nvarcalc()<<endl;
    R_report   << endl;
    
  R_report<<"$Steep_Prior" <<endl;
  R_report<<steepnessprior<<" "<<
    cvsteepnessprior<<" "<<
    phase_srec<<" "<< endl;
    R_report   << endl;

  R_report<<"$sigmarPrior " <<endl;
  R_report<<sigmarprior<<" "<<  cvsigmarprior <<" "<<phase_sigmar<<endl;
  R_report   << endl;

  R_report<<"$Rec_estimated_in_styr_endyr " <<endl;
  R_report<<styr_rec    <<" "<<endyr        <<" "<<endl;
  R_report   << endl;
  R_report<<"$SR_Curve_fit__in_styr_endyr " <<endl;
  R_report<<styr_rec_est<<" "<<endyr_rec_est<<" "<<endl;
  R_report   << endl;
  R_report<<"$Model_styr_endyr" <<endl;
  R_report<<styr        <<" "<<endyr        <<" "<<endl;
  R_report   << endl;

  R_report<<"$M_prior "<<endl;
  R_report<< natmortprior<< " "<< cvnatmortprior<<" "<<phase_M<<endl;
  R_report   << endl;
  R_report<<"$qprior " <<endl;
  R_report<< qprior<<" "<<cvqprior<<" "<< phase_q<<endl;
  R_report<<"$q_power_prior " <<endl;
  R_report<< q_power_prior<<" "<<cvq_power_prior<<" "<< phase_q_power<<endl;
  R_report   << endl;

  R_report<<"$cv_catchbiomass " <<endl;
  R_report<<cv_catchbiomass<<" "<<endl;
  R_report   << endl;
  R_report<<"$Projection_years"<<endl;
  R_report<< nproj_yrs<<endl;
  R_report   << endl;
  
  R_report << "$Fsh_sel_opt_fish "<<endl;
  for (int k=1;k<=nfsh;k++)
    R_report<<k<<" "<<fsh_sel_opt(k)<<" "<<sel_change_in_fsh(k)<<endl;
    R_report   << endl;
   R_report<<"$Index_Sel_Opt" <<endl;
  for (int k=1;k<=nind;k++)
  R_report<<k<<" "<<(ind_sel_opt(k))<<endl;
  R_report   << endl;
    
  R_report <<"$Phase_index_Sel_Coffs "<<endl;
  R_report <<phase_selcoff_ind<<endl;
  R_report   << endl;
  R_report <<"$Fshry_Selages " << endl;
  R_report << nselages_in_fsh  <<endl;
  R_report   << endl;
  R_report <<"$Index_Selages " <<endl;
  R_report <<nselages_in_ind <<endl;
  R_report   << endl;

  R_report << "$Phase_for_age_spec_fishery"<<endl;
  R_report <<phase_selcoff_fsh<<endl;
  R_report   << endl;
  R_report << "$Phase_for_logistic_fishery"<<endl;
  R_report <<phase_logist_fsh<<endl;
  R_report   << endl;
  R_report << "$Phase_for_dble_logistic_fishery "<<endl;
  R_report <<phase_dlogist_fsh<<endl;
  R_report   << endl;

  R_report << "$Phase_for_age_spec_index  "<<endl;
  R_report <<phase_selcoff_ind<<endl;
  R_report   << endl;
  R_report << "$Phase_for_logistic_index  "<<endl;
  R_report <<phase_logist_ind<<endl;
  R_report   << endl;
  R_report << "$Phase_for_dble_logistic_indy "<<endl;
  R_report <<phase_dlogist_ind<<endl;
  R_report   << endl;
  
  for (int k=1;k<=nfsh;k++)
  {
    if (nyrs_fsh_age(k)>0)
    {
      R_report <<"$EffN_Fsh_"<<(k)<<""<<endl;
      for (i=1;i<=nyrs_fsh_age(k);i++)
      {
        double sda_tmp = Sd_age(oac_fsh(k,i));
        R_report << yrs_fsh_age(k,i);
        R_report<< " "<<Eff_N(oac_fsh(k,i),eac_fsh(k,i)) ;
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

  for (int k=1;k<=nfsh;k++)
  {
    R_report <<"$C_fsh_" <<(k)<<"" << endl; 
    for (i=styr;i<=endyr;i++)
      R_report <<i<<" "<<catage(k,i)<< endl;
  }

  R_report <<"$wt_a_pop" << endl<< wt_pop  <<endl;
  R_report <<"$mature_a" << endl<< maturity<<endl;
  for (int k=1;k<=nfsh;k++)
  {
    R_report <<"$wt_fsh_"<<(k)<<""<<endl;
    for (i=styr;i<=endyr;i++)
      R_report <<i<<" "<<wt_fsh(k,i)<< endl;
  }
  
  for (int k=1;k<=nind;k++)
  {
    R_report <<"$wt_ind_"<<(k)<<""<<endl;
    for (i=styr;i<=endyr;i++)
      R_report <<i<<" "<<wt_ind(k,i)<< endl;
  }
  for (int k=1;k<=nind;k++)
  {
    if (nyrs_ind_age(k)>0)
    {
      R_report <<"$EffN_Index_"<<(k)<<""<<endl;
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

  R_report.close();


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

FUNCTION double get_AC(const int& indind)
  // Functions to compute autocorrelation in residuals 
  int i1,i2,yy;
  i1 = 1;
  i2 = nyrs_ind(indind);
  double actmp;
  dvector res(1,i2);
  for (i=1;i<=i2;i++)
  {
    yy = int(yrs_ind(indind,i));
    cout<<yy<<" "<<obs_ind(indind,i)<<" " <<pred_ind(indind,yy)<<endl;
    res(i) = log(obs_ind(indind,i)) - value(log(pred_ind(indind,yy)));
  }
  double m1 = (mean(res(i1,i2-1)));
  double m2 = (mean(res(i1+1,i2))); 
  actmp = mean( elem_prod( ++res(i1,i2-1) - m1, res(i1+1,i2) - m2)) /
          (sqrt(mean( square(res(i1,i2-1) - m1 )))  * sqrt(mean(square(res(i1+1,i2) - m2 ))) );
  return(actmp);
