#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  extern "C"
  {
    void do_not_exit(const int exit_code)
    {
       // Does Nothing
    }
  }
  #include <float.h>
  #include <admodel.h>
  #undef COUT
  #define COUT(object) cout << #object "," << object << endl;
  ofstream write_log("Input_Log.rep");
  #undef write_log
  #define write_log(object) write_log << #object "\n" << object << endl;
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <gp.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
 int on,opt;
 retro_yrs=0;
 isretro=0;
 if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1){
   retro_yrs=atoi(ad_comm::argv[on+1]);
   isretro=1;
   cout << "\n!!  Implementing retro w/ peels="<< retro_yrs << " !!\n";
 }
  pad_report1 = new ofstream("mceval.dat");
  styr.allocate("styr");
  endyr.allocate("endyr");
  rcrage.allocate("rcrage");
  trmage.allocate("trmage");
  nbins1.allocate("nbins1");
  nbins2.allocate("nbins2");
  nbins3.allocate("nbins3");
 write_log(nbins3);
  cattot.allocate(styr,endyr,"cattot");
  cattot_log_sd.allocate(styr,endyr,"cattot_log_sd");
 write_log(cattot_log_sd);
  nyrs_fsh.allocate("nyrs_fsh");
  fshyrs.allocate(1,nyrs_fsh,"fshyrs");
  multN_fsh.allocate(1,nyrs_fsh,"multN_fsh");
  ac_yng_fsh.allocate(1,nyrs_fsh,"ac_yng_fsh");
  ac_old_fsh.allocate(1,nyrs_fsh,"ac_old_fsh");
  nyrslen_fsh.allocate("nyrslen_fsh");
 write_log(nyrslen_fsh);
  fshlenyrs.allocate(1,nyrslen_fsh,"fshlenyrs");
 write_log(fshlenyrs);
  multNlen_fsh.allocate(1,nyrslen_fsh,"multNlen_fsh");
 write_log(multNlen_fsh);
  rwlk_sd.allocate(styr,endyr-1,"rwlk_sd");
  catp.allocate(1,nyrs_fsh,rcrage,trmage,"catp");
 write_log(catp);
  lenp.allocate(1,nyrslen_fsh,1,nbins1,"lenp");
 write_log(lenp);
  wt_fsh.allocate(styr,endyr,rcrage,trmage,"wt_fsh");
 write_log(wt_fsh);
  nyrs_srv1.allocate("nyrs_srv1");
 write_log(nyrs_srv1);
  srvyrs1.allocate(1,nyrs_srv1,"srvyrs1");
  indxsurv1.allocate(1,nyrs_srv1,"indxsurv1");
  indxsurv_log_sd1.allocate(1,nyrs_srv1,"indxsurv_log_sd1");
 write_log(indxsurv_log_sd1);
  q1_rwlk_sd.allocate(styr,endyr-1,"q1_rwlk_sd");
 write_log(q1_rwlk_sd);
  yrfrct_srv1.allocate(styr,endyr,"yrfrct_srv1");
 write_log(yrfrct_srv1);
  nyrsac_srv1.allocate("nyrsac_srv1");
  srv_acyrs1.allocate(1,nyrsac_srv1,"srv_acyrs1");
  multN_srv1.allocate(1,nyrsac_srv1,"multN_srv1");
 write_log(multN_srv1);
  ac_yng_srv1.allocate(1,nyrsac_srv1,"ac_yng_srv1");
  ac_old_srv1.allocate(1,nyrsac_srv1,"ac_old_srv1");
  nyrslen_srv1.allocate("nyrslen_srv1");
 write_log(nyrslen_srv1);
  srv_lenyrs1.allocate(1,nyrslen_srv1,"srv_lenyrs1");
  multNlen_srv1.allocate(1,nyrslen_srv1,"multNlen_srv1");
 write_log(multNlen_srv1);
  srvp1.allocate(1,nyrsac_srv1,rcrage,trmage,"srvp1");
 write_log(srvp1);
  srvlenp1.allocate(1,nyrslen_srv1,1,nbins3,"srvlenp1");
 write_log(srvlenp1);
  wt_srv1.allocate(styr,endyr,rcrage,trmage,"wt_srv1");
 write_log(wt_srv1);
  nyrs_srv2.allocate("nyrs_srv2");
 write_log(nyrs_srv2);
  srvyrs2.allocate(1,nyrs_srv2,"srvyrs2");
 write_log(srvyrs2);
  indxsurv2.allocate(1,nyrs_srv2,"indxsurv2");
  indxsurv_log_sd2.allocate(1,nyrs_srv2,"indxsurv_log_sd2");
 write_log(indxsurv2);
  q2_rwlk_sd.allocate(styr,endyr-1,"q2_rwlk_sd");
 write_log(q2_rwlk_sd);
  yrfrct_srv2.allocate(styr,endyr,"yrfrct_srv2");
  nyrsac_srv2.allocate("nyrsac_srv2");
  srv_acyrs2.allocate(1,nyrsac_srv2,"srv_acyrs2");
  multN_srv2.allocate(1,nyrsac_srv2,"multN_srv2");
 write_log(multN_srv2);
  ac_yng_srv2.allocate(1,nyrsac_srv2,"ac_yng_srv2");
  ac_old_srv2.allocate(1,nyrsac_srv2,"ac_old_srv2");
  nyrslen_srv2.allocate("nyrslen_srv2");
  srv_lenyrs2.allocate(1,nyrslen_srv2,"srv_lenyrs2");
  multNlen_srv2.allocate(1,nyrslen_srv2,"multNlen_srv2");
  srvp2.allocate(1,nyrsac_srv2,rcrage,trmage,"srvp2");
  srvlenp2.allocate(1,nyrslen_srv2,1,nbins2,"srvlenp2");
  wt_srv2.allocate(styr,endyr,rcrage,trmage,"wt_srv2");
 write_log(wt_srv2);
  nyrs_srv3.allocate("nyrs_srv3");
  srvyrs3.allocate(1,nyrs_srv3,"srvyrs3");
  indxsurv3.allocate(1,nyrs_srv3,"indxsurv3");
  indxsurv_log_sd3.allocate(1,nyrs_srv3,"indxsurv_log_sd3");
  q3_rwlk_sd.allocate(styr,endyr-1,"q3_rwlk_sd");
 write_log(indxsurv3);
  yrfrct_srv3.allocate(styr,endyr,"yrfrct_srv3");
  write_log(yrfrct_srv3);
  nyrsac_srv3.allocate("nyrsac_srv3");
 write_log(nyrsac_srv3);
  srv_acyrs3.allocate(1,nyrsac_srv3,"srv_acyrs3");
  write_log(srv_acyrs3);
  multN_srv3.allocate(1,nyrsac_srv3,"multN_srv3");
 write_log(multN_srv3);
  nyrslen_srv3.allocate("nyrslen_srv3");
  srv_lenyrs3.allocate(1,nyrslen_srv3,"srv_lenyrs3");
  multNlen_srv3.allocate(1,nyrslen_srv3,"multNlen_srv3");
  srvp3.allocate(1,nyrsac_srv3,rcrage,trmage,"srvp3");
 write_log(srvp3);
  srvlenp3.allocate(1,nyrslen_srv3,1,nbins2,"srvlenp3");
 write_log(srvlenp3);
  wt_srv3.allocate(styr,endyr,rcrage,trmage,"wt_srv3");
 write_log(wt_srv3);
  nyrs_srv4.allocate("nyrs_srv4");
 write_log(nyrs_srv4);
  srvyrs4.allocate(1,nyrs_srv4,"srvyrs4");
  indxsurv4.allocate(1,nyrs_srv4,"indxsurv4");
 write_log(indxsurv4);
  indxsurv_log_sd4.allocate(1,nyrs_srv4,"indxsurv_log_sd4");
  nyrs_srv5.allocate("nyrs_srv5");
 write_log(nyrs_srv5);
  srvyrs5.allocate(1,nyrs_srv5,"srvyrs5");
  indxsurv5.allocate(1,nyrs_srv5,"indxsurv5");
  indxsurv_log_sd5.allocate(1,nyrs_srv5,"indxsurv_log_sd5");
  nyrs_srv6.allocate("nyrs_srv6");
 write_log(nyrs_srv6);
  srvyrs6.allocate(1,nyrs_srv6,"srvyrs6");
 write_log(srvyrs6);
  indxsurv6.allocate(1,nyrs_srv6,"indxsurv6");
 write_log(nyrsac_srv6);
  indxsurv_log_sd6.allocate(1,nyrs_srv6,"indxsurv_log_sd6");
  yrfrct_srv6.allocate(styr,endyr,"yrfrct_srv6");
 write_log(yrfrct_srv6);
  nyrsac_srv6.allocate("nyrsac_srv6");
 write_log(nyrsac_srv6);
  srv_acyrs6.allocate(1,nyrsac_srv6,"srv_acyrs6");
  multN_srv6.allocate(1,nyrsac_srv6,"multN_srv6");
 write_log(multN_srv6);
  ac_yng_srv6.allocate(1,nyrsac_srv6,"ac_yng_srv6");
 write_log(ac_yng_srv6);
  ac_old_srv6.allocate(1,nyrsac_srv6,"ac_old_srv6");
  nyrslen_srv6.allocate("nyrslen_srv6");
 write_log(nyrslen_srv6);
  srv_lenyrs6.allocate(1,nyrslen_srv6,"srv_lenyrs6");
  multNlen_srv6.allocate(1,nyrslen_srv6,"multNlen_srv6");
 write_log(multNlen_srv6);
  srvp6.allocate(1,nyrsac_srv6,rcrage,trmage,"srvp6");
 write_log(srvp6);
  srvlenp6.allocate(1,nyrslen_srv6,1,nbins2,"srvlenp6");
  wt_srv6.allocate(styr,endyr,rcrage,trmage,"wt_srv6");
 write_log(wt_srv6);
  age_trans.allocate(rcrage,trmage,rcrage,trmage,"age_trans");
 write_log(age_trans);
  len_trans1.allocate(rcrage,trmage,1,nbins1,"len_trans1");
 write_log(len_trans1);
  len_trans2.allocate(rcrage,trmage,1,nbins2,"len_trans2");
  len_trans3.allocate(rcrage,trmage,1,nbins3,"len_trans3");
  wt_pop.allocate(styr,endyr,rcrage,trmage);
  wt_spawn.allocate(styr,endyr,rcrage,trmage);
 wt_pop=wt_srv2;
 wt_spawn=wt_srv1;
  mat.allocate(rcrage,trmage,"mat");
 write_log(mat);
  wt_pop_proj.allocate(rcrage,trmage);
  wt_spawn_proj.allocate(rcrage,trmage);
  wt_fsh_proj.allocate(rcrage,trmage);
  wt_srv_proj.allocate(rcrage,trmage);
   // for projections take averages of WAA only from recent survey years with data
    wt_pop_proj.initialize();
    wt_spawn_proj.initialize();
    for(int a=rcrage;a<=trmage;a++){
      for(int i=1;i<=3;i++)
          wt_pop_proj(a)+=wt_srv2(srv_acyrs2(nyrsac_srv2-i+1),a)/3;
      for(int i=1;i<=5;i++){
          wt_spawn_proj(a)+=wt_srv1(srv_acyrs1(nyrsac_srv1-i+1),a)/5;
      // for predicting what the survey will see next year, Shelikof for now
          wt_srv_proj(a)+=wt_srv1(srv_acyrs1(nyrsac_srv1-i+1),a)/5;
      }	  
     }
     wt_fsh_proj=wt_fsh(endyr);
  Ftarget.allocate(endyr+1,endyr+5,"Ftarget");
  B40.allocate("B40");
  log_mean_recr_proj.allocate("log_mean_recr_proj");
  sigmasq_recr.allocate("sigmasq_recr");
  check.allocate("check");
 if(check != -999){ cerr << "Failed to parse .dat file, final value=" << check <<endl; ad_exit(1);}
  if(retro_yrs<0){cerr << "bad peels in -retro option" << endl; ad_exit(1);};
  endyr0=endyr;
  endyr=endyr-retro_yrs;
}

void model_parameters::initializationfunction(void)
{
  mean_log_initN.set_initial_value(0.0);
  mean_log_recruit.set_initial_value(12.0);
  mean_log_F.set_initial_value(-1.6);
  M.set_initial_value(0.30);
  log_q1_mean.set_initial_value(0.0);
  log_q1_dev.set_initial_value(0.0);
  log_q2_mean.set_initial_value(0.0);
  log_q2_dev.set_initial_value(0.0);
  log_q3_mean.set_initial_value(-1.6);
  log_q3_dev.set_initial_value(0.0);
  log_q6.set_initial_value(0.0);
  slp1_fsh_dev.set_initial_value(0.0);
  inf1_fsh_dev.set_initial_value(0.0);
  slp2_fsh_dev.set_initial_value(0.0);
  inf2_fsh_dev.set_initial_value(0.0);
  log_slp1_fsh_mean.set_initial_value(1.0);
  inf1_fsh_mean.set_initial_value(4.0);
  log_slp2_fsh_mean.set_initial_value(1.0);
  inf2_fsh_mean.set_initial_value(8.0);
  log_slp2_srv1.set_initial_value(0.5);
  inf2_srv1.set_initial_value(9.0);
  log_slp1_srv2.set_initial_value(-0.8);
  inf1_srv2.set_initial_value(4.0);
  log_slp2_srv2.set_initial_value(1);
  inf2_srv2.set_initial_value(20);
  log_slp1_srv3.set_initial_value(0.0);
  inf1_srv3.set_initial_value(5.0);
  log_slp1_srv6.set_initial_value(4.9);
  inf1_srv6.set_initial_value(0.5);
  log_slp2_srv6.set_initial_value(1);
  inf2_srv6.set_initial_value(20);
  sigmaR.set_initial_value(1.3);
  natMscalar.set_initial_value(1);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  M.allocate(rcrage,trmage,0.1,5.0,-1,"M");
  mean_log_initN.allocate(-15,15,-1,"mean_log_initN");
  dev_log_initN.allocate(rcrage+1,trmage,-15,15,-2,"dev_log_initN");
  initN.allocate(rcrage+1,trmage,"initN");
  #ifndef NO_AD_INITIALIZE
    initN.initialize();
  #endif
  mean_log_recruit.allocate(-15,15,1,"mean_log_recruit");
  dev_log_recruit.allocate(styr,endyr,-15,15,3,"dev_log_recruit");
  sigmaR.allocate(0,5,-1,"sigmaR");
  recruit.allocate(styr,endyr,"recruit");
  log_recruit.allocate(styr,endyr,"log_recruit");
  log_recr_proj.allocate(endyr+1,endyr+5,-5,5,10,"log_recr_proj");
  recruit_proj.allocate(endyr+1,endyr+5,"recruit_proj");
  N_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  F_proj.allocate(endyr+1,endyr+5,"F_proj");
  #ifndef NO_AD_INITIALIZE
    F_proj.initialize();
  #endif
  Z_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  C_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"C_proj");
  #ifndef NO_AD_INITIALIZE
    C_proj.initialize();
  #endif
  Nsrv_proj.allocate(endyr+1,endyr+5,rcrage,trmage,"Nsrv_proj");
  #ifndef NO_AD_INITIALIZE
    Nsrv_proj.initialize();
  #endif
  slctfsh_proj.allocate(rcrage,trmage,"slctfsh_proj");
  #ifndef NO_AD_INITIALIZE
    slctfsh_proj.initialize();
  #endif
  Ecattot_proj.allocate(endyr+1,endyr+5,"Ecattot_proj");
  #ifndef NO_AD_INITIALIZE
    Ecattot_proj.initialize();
  #endif
  Esumbio_proj.allocate(endyr+1,endyr+5,"Esumbio_proj");
  Espawnbio_proj.allocate(endyr+1,endyr+5,"Espawnbio_proj");
  Esrv_proj.allocate(endyr+1,endyr+5,"Esrv_proj");
  Exrate_proj.allocate(endyr+1,endyr+5,"Exrate_proj");
  sbio.allocate("sbio");
  #ifndef NO_AD_INITIALIZE
  sbio.initialize();
  #endif
  log_slp1_fsh_mean.allocate(-5,5,4,"log_slp1_fsh_mean");
  inf1_fsh_mean.allocate(1,5,4,"inf1_fsh_mean");
  log_slp2_fsh_mean.allocate(-5,5,4,"log_slp2_fsh_mean");
  inf2_fsh_mean.allocate(7,20,6,"inf2_fsh_mean");
  slp1_fsh_dev.allocate(styr,endyr,-5,5,7,"slp1_fsh_dev");
  inf1_fsh_dev.allocate(styr,endyr,-5,5,7,"inf1_fsh_dev");
  slp2_fsh_dev.allocate(styr,endyr,-5,5,6,"slp2_fsh_dev");
  inf2_fsh_dev.allocate(styr,endyr,-5,5,7,"inf2_fsh_dev");
  slp1_fsh.allocate(styr,endyr,"slp1_fsh");
  #ifndef NO_AD_INITIALIZE
    slp1_fsh.initialize();
  #endif
  inf1_fsh.allocate(styr,endyr,"inf1_fsh");
  #ifndef NO_AD_INITIALIZE
    inf1_fsh.initialize();
  #endif
  slp2_fsh.allocate(styr,endyr,"slp2_fsh");
  #ifndef NO_AD_INITIALIZE
    slp2_fsh.initialize();
  #endif
  inf2_fsh.allocate(styr,endyr,"inf2_fsh");
  #ifndef NO_AD_INITIALIZE
    inf2_fsh.initialize();
  #endif
  log_slp2_srv1.allocate(-5,5,7,"log_slp2_srv1");
  inf2_srv1.allocate(3,10,7,"inf2_srv1");
  log_slp1_srv2.allocate(-5,5,6,"log_slp1_srv2");
  inf1_srv2.allocate(1,50,6,"inf1_srv2");
  log_slp2_srv2.allocate(-5,5,-1,"log_slp2_srv2");
  inf2_srv2.allocate(5,25,-1,"inf2_srv2");
  log_slp1_srv3.allocate(-5,5,9,"log_slp1_srv3");
  inf1_srv3.allocate(1,20,9,"inf1_srv3");
  log_slp1_srv6.allocate(-5,5,-1,"log_slp1_srv6");
  inf1_srv6.allocate(0,50,-1,"inf1_srv6");
  log_slp2_srv6.allocate(-5,5,7,"log_slp2_srv6");
  inf2_srv6.allocate(5,25,7,"inf2_srv6");
  mean_log_F.allocate(-10,10,1,"mean_log_F");
  dev_log_F.allocate(styr,endyr,-10,10,2,"dev_log_F");
  F.allocate(styr,endyr,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  log_q1_mean.allocate(-10,10,5,"log_q1_mean");
  log_q1_dev.allocate(styr,endyr,-5,5,-5,"log_q1_dev");
  log_q2_mean.allocate(-10,10,5,"log_q2_mean");
  log_q2_dev.allocate(styr,endyr,-5,5,-1,"log_q2_dev");
  log_q3_mean.allocate(-10,10,6,"log_q3_mean");
  log_q3_dev.allocate(styr,endyr,-5,5,5,"log_q3_dev");
  log_q4.allocate(-10,10,6,"log_q4");
  q4_pow.allocate(-10,10,-6,"q4_pow");
  log_q5.allocate(-10,10,6,"log_q5");
  q5_pow.allocate(-10,10,-6,"q5_pow");
  log_q6.allocate(-10,10,5,"log_q6");
  natMscalar.allocate(0,5,-5,"natMscalar");
  q1.allocate(styr,endyr,"q1");
  #ifndef NO_AD_INITIALIZE
    q1.initialize();
  #endif
  q2.allocate(styr,endyr,"q2");
  #ifndef NO_AD_INITIALIZE
    q2.initialize();
  #endif
  q3.allocate(styr,endyr,"q3");
  #ifndef NO_AD_INITIALIZE
    q3.initialize();
  #endif
  log_q1.allocate(styr,endyr,"log_q1");
  log_q2.allocate(styr,endyr,"log_q2");
  log_q3.allocate(styr,endyr,"log_q3");
  q4.allocate("q4");
  #ifndef NO_AD_INITIALIZE
  q4.initialize();
  #endif
  q5.allocate("q5");
  #ifndef NO_AD_INITIALIZE
  q5.initialize();
  #endif
  q6.allocate("q6");
  #ifndef NO_AD_INITIALIZE
  q6.initialize();
  #endif
  N.allocate(styr,endyr,rcrage,trmage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  endN.allocate(rcrage,trmage,"endN");
  Z.allocate(styr,endyr,rcrage,trmage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  C.allocate(styr,endyr,rcrage,trmage,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  Nsrv1.allocate(styr,endyr,rcrage,trmage,"Nsrv1");
  #ifndef NO_AD_INITIALIZE
    Nsrv1.initialize();
  #endif
  slctsrv1.allocate(rcrage,trmage,"slctsrv1");
  Nsrv2.allocate(styr,endyr,rcrage,trmage,"Nsrv2");
  #ifndef NO_AD_INITIALIZE
    Nsrv2.initialize();
  #endif
  slctsrv2.allocate(rcrage,trmage,"slctsrv2");
  Nsrv3.allocate(styr,endyr,rcrage,trmage,"Nsrv3");
  #ifndef NO_AD_INITIALIZE
    Nsrv3.initialize();
  #endif
  slctsrv3.allocate(rcrage,trmage,"slctsrv3");
  Nsrv6.allocate(styr,endyr,rcrage,trmage,"Nsrv6");
  #ifndef NO_AD_INITIALIZE
    Nsrv6.initialize();
  #endif
  slctsrv6.allocate(rcrage,trmage,"slctsrv6");
  #ifndef NO_AD_INITIALIZE
    slctsrv6.initialize();
  #endif
  slctfsh.allocate(styr,endyr,rcrage,trmage,"slctfsh");
  #ifndef NO_AD_INITIALIZE
    slctfsh.initialize();
  #endif
  slctsrv1_logit.allocate(rcrage,trmage,"slctsrv1_logit");
  slctsrv2_logit.allocate(rcrage,trmage,"slctsrv2_logit");
  slctsrv3_logit.allocate(rcrage,trmage,"slctsrv3_logit");
  slctsrv6_logit.allocate(rcrage,trmage,"slctsrv6_logit");
  slctfsh_logit.allocate(rcrage,trmage,"slctfsh_logit");
  Eecocon.allocate(styr,endyr,"Eecocon");
  #ifndef NO_AD_INITIALIZE
    Eecocon.initialize();
  #endif
  Eec.allocate(styr,endyr,rcrage,trmage,"Eec");
  #ifndef NO_AD_INITIALIZE
    Eec.initialize();
  #endif
  Ecattot.allocate(styr,endyr,"Ecattot");
  #ifndef NO_AD_INITIALIZE
    Ecattot.initialize();
  #endif
  Ecatp.allocate(styr,endyr,rcrage,trmage,"Ecatp");
  #ifndef NO_AD_INITIALIZE
    Ecatp.initialize();
  #endif
  Elenp.allocate(styr,endyr,1,nbins1,"Elenp");
  #ifndef NO_AD_INITIALIZE
    Elenp.initialize();
  #endif
  Eindxsurv1.allocate(styr,endyr,"Eindxsurv1");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv1.initialize();
  #endif
  Esrvp1.allocate(styr,endyr,rcrage,trmage,"Esrvp1");
  #ifndef NO_AD_INITIALIZE
    Esrvp1.initialize();
  #endif
  Esrvlenp1.allocate(styr,endyr,1,nbins3,"Esrvlenp1");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp1.initialize();
  #endif
  Eindxsurv2.allocate(styr,endyr,"Eindxsurv2");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv2.initialize();
  #endif
  Esrvp2.allocate(styr,endyr,rcrage,trmage,"Esrvp2");
  #ifndef NO_AD_INITIALIZE
    Esrvp2.initialize();
  #endif
  Esrvlenp2.allocate(styr,endyr,1,nbins2,"Esrvlenp2");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp2.initialize();
  #endif
  Eindxsurv3.allocate(styr,endyr,"Eindxsurv3");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv3.initialize();
  #endif
  Esrvp3.allocate(styr,endyr,rcrage,trmage,"Esrvp3");
  #ifndef NO_AD_INITIALIZE
    Esrvp3.initialize();
  #endif
  Esrvlenp3.allocate(styr,endyr,1,nbins2,"Esrvlenp3");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp3.initialize();
  #endif
  Eindxsurv4.allocate(styr,endyr,"Eindxsurv4");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv4.initialize();
  #endif
  Eindxsurv5.allocate(styr,endyr,"Eindxsurv5");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv5.initialize();
  #endif
  Eindxsurv6.allocate(styr,endyr,"Eindxsurv6");
  #ifndef NO_AD_INITIALIZE
    Eindxsurv6.initialize();
  #endif
  Esrvp6.allocate(styr,endyr,rcrage,trmage,"Esrvp6");
  #ifndef NO_AD_INITIALIZE
    Esrvp6.initialize();
  #endif
  Esrvlenp6.allocate(styr,endyr,1,nbins2,"Esrvlenp6");
  #ifndef NO_AD_INITIALIZE
    Esrvlenp6.initialize();
  #endif
  Espawnbio.allocate(styr,endyr,"Espawnbio");
  Esumbio.allocate(styr,endyr,"Esumbio");
  Espawnbio_log.allocate(styr,endyr,"Espawnbio_log");
  Esumbio_log.allocate(styr,endyr,"Esumbio_log");
  Espawnbio_2plus.allocate(styr,endyr,"Espawnbio_2plus");
  #ifndef NO_AD_INITIALIZE
    Espawnbio_2plus.initialize();
  #endif
  Etotalbio.allocate(styr,endyr,"Etotalbio");
  #ifndef NO_AD_INITIALIZE
    Etotalbio.initialize();
  #endif
  loglik.allocate(1,24,"loglik");
  #ifndef NO_AD_INITIALIZE
    loglik.initialize();
  #endif
  llcatp.allocate(1,nyrs_fsh,"llcatp");
  #ifndef NO_AD_INITIALIZE
    llcatp.initialize();
  #endif
  lllenp.allocate(1,nyrslen_fsh,"lllenp");
  #ifndef NO_AD_INITIALIZE
    lllenp.initialize();
  #endif
  llsrvp1.allocate(1,nyrsac_srv1,"llsrvp1");
  #ifndef NO_AD_INITIALIZE
    llsrvp1.initialize();
  #endif
  llsrvlenp1.allocate(1,nyrslen_srv1,"llsrvlenp1");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp1.initialize();
  #endif
  llsrvp2.allocate(1,nyrsac_srv2,"llsrvp2");
  #ifndef NO_AD_INITIALIZE
    llsrvp2.initialize();
  #endif
  llsrvlenp2.allocate(1,nyrslen_srv2,"llsrvlenp2");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp2.initialize();
  #endif
  llsrvp3.allocate(1,nyrsac_srv3,"llsrvp3");
  #ifndef NO_AD_INITIALIZE
    llsrvp3.initialize();
  #endif
  llsrvlenp3.allocate(1,nyrslen_srv3,"llsrvlenp3");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp3.initialize();
  #endif
  llsrvp6.allocate(1,nyrsac_srv6,"llsrvp6");
  #ifndef NO_AD_INITIALIZE
    llsrvp6.initialize();
  #endif
  llsrvlenp6.allocate(1,nyrslen_srv6,"llsrvlenp6");
  #ifndef NO_AD_INITIALIZE
    llsrvlenp6.initialize();
  #endif
  res_fish.allocate(1,nyrs_fsh,rcrage,2*trmage-rcrage+1,"res_fish");
  #ifndef NO_AD_INITIALIZE
    res_fish.initialize();
  #endif
  res_srv1.allocate(1,nyrsac_srv1,rcrage,2*trmage-rcrage+1,"res_srv1");
  #ifndef NO_AD_INITIALIZE
    res_srv1.initialize();
  #endif
  res_srv2.allocate(1,nyrsac_srv2,rcrage,2*trmage-rcrage+1,"res_srv2");
  #ifndef NO_AD_INITIALIZE
    res_srv2.initialize();
  #endif
  res_srv3.allocate(1,nyrsac_srv3,rcrage,2*trmage-rcrage+1,"res_srv3");
  #ifndef NO_AD_INITIALIZE
    res_srv3.initialize();
  #endif
  res_srv3len.allocate(1,nyrslen_srv3,1,2*nbins2,"res_srv3len");
  #ifndef NO_AD_INITIALIZE
    res_srv3len.initialize();
  #endif
  res_srv6.allocate(1,nyrsac_srv6,rcrage,2*trmage-rcrage+1,"res_srv6");
  #ifndef NO_AD_INITIALIZE
    res_srv6.initialize();
  #endif
  pearson_fish.allocate(1,nyrs_fsh,rcrage,trmage,"pearson_fish");
  #ifndef NO_AD_INITIALIZE
    pearson_fish.initialize();
  #endif
  pearson_srv1.allocate(1,nyrsac_srv1,rcrage,trmage,"pearson_srv1");
  #ifndef NO_AD_INITIALIZE
    pearson_srv1.initialize();
  #endif
  pearson_srv2.allocate(1,nyrsac_srv2,rcrage,trmage,"pearson_srv2");
  #ifndef NO_AD_INITIALIZE
    pearson_srv2.initialize();
  #endif
  pearson_srv3.allocate(1,nyrsac_srv3,rcrage,trmage,"pearson_srv3");
  #ifndef NO_AD_INITIALIZE
    pearson_srv3.initialize();
  #endif
  pearson_srv3len.allocate(1,nyrslen_srv3,1,nbins2,"pearson_srv3len");
  #ifndef NO_AD_INITIALIZE
    pearson_srv3len.initialize();
  #endif
  pearson_srv6.allocate(1,nyrsac_srv6,rcrage,trmage,"pearson_srv6");
  #ifndef NO_AD_INITIALIZE
    pearson_srv6.initialize();
  #endif
  effN_fsh.allocate(1,nyrs_fsh,"effN_fsh");
  #ifndef NO_AD_INITIALIZE
    effN_fsh.initialize();
  #endif
  effN_srv1.allocate(1,nyrsac_srv1,"effN_srv1");
  #ifndef NO_AD_INITIALIZE
    effN_srv1.initialize();
  #endif
  effN_srv2.allocate(1,nyrsac_srv2,"effN_srv2");
  #ifndef NO_AD_INITIALIZE
    effN_srv2.initialize();
  #endif
  effN_srv3.allocate(1,nyrsac_srv3,"effN_srv3");
  #ifndef NO_AD_INITIALIZE
    effN_srv3.initialize();
  #endif
  effN_srv6.allocate(1,nyrsac_srv6,"effN_srv6");
  #ifndef NO_AD_INITIALIZE
    effN_srv6.initialize();
  #endif
  RMSE_srv1.allocate("RMSE_srv1");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv1.initialize();
  #endif
  RMSE_srv2.allocate("RMSE_srv2");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv2.initialize();
  #endif
  RMSE_srv3.allocate("RMSE_srv3");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv3.initialize();
  #endif
  RMSE_srv4.allocate("RMSE_srv4");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv4.initialize();
  #endif
  RMSE_srv5.allocate("RMSE_srv5");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv5.initialize();
  #endif
  RMSE_srv6.allocate("RMSE_srv6");
  #ifndef NO_AD_INITIALIZE
  RMSE_srv6.initialize();
  #endif
  var_prof.allocate("var_prof");
  objfun.allocate("objfun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 // Do upper and lower accumulations in age composition data
 // Fishery
  for (i=1;i<=nyrs_fsh;i++) {
   for (j=rcrage;j<=trmage;j++) {
     if(j<ac_yng_fsh(i)) {
       catp(i,ac_yng_fsh(i)) += catp(i,j);
       catp(i,j) = 0;
     }
     if(j>ac_old_fsh(i)) {
       catp(i,ac_old_fsh(i)) += catp(i,j);
       catp(i,j) = 0;
     }
    }
		catp(i) = catp(i)/sum(catp(i));
  }
 for (i=1;i<=nyrsac_srv1;i++) {
   for (j=rcrage;j<=trmage;j++) {
     if(j<ac_yng_srv1(i)) {
       srvp1(i,ac_yng_srv1(i)) += srvp1(i,j);
       srvp1(i,j) = 0;
     }
     if(j>ac_old_srv1(i)) {
       srvp1(i,ac_old_srv1(i)) += srvp1(i,j);
       srvp1(i,j) = 0;
     }
   }
		srvp1(i) = srvp1(i)/sum(srvp1(i));
 }
  for (i=1;i<=nyrsac_srv2;i++) {
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv2(i)) {
	srvp2(i,ac_yng_srv2(i)) += srvp2(i,j);
	srvp2(i,j) = 0;
      }
      if(j>ac_old_srv2(i)) {
	srvp2(i,ac_old_srv2(i)) += srvp2(i,j);
	srvp2(i,j) = 0;
      }
    }
		srvp2(i) = srvp2(i)/sum(srvp2(i));
  }
	
  // Survey 6
  for (i=1;i<=nyrsac_srv6;i++) {
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv6(i)) {
	srvp6(i,ac_yng_srv6(i)) += srvp6(i,j);
	srvp6(i,j) = 0;
      }
      if(j>ac_old_srv6(i)) {
	srvp6(i,ac_old_srv6(i)) += srvp6(i,j);
	srvp6(i,j) = 0;
      }
    }
		srvp6(i) = srvp6(i)/sum(srvp6(i));
  }
  o = 0.00001;
  var_prof.set_stepnumber(30);
  //  var_prof.set_stepnumber(10);
  var_prof.set_stepsize(0.1);
}

void model_parameters::userfunction(void)
{
  objfun =0.0;
  ofstream& report1= *pad_report1;
  Convert_log_parameters();
  Selectivity();
  Mortality();
  Numbers_at_age();
  Catch_at_age();
  Expected_values();
  if(last_phase())
  {
    Projections();
  } 
  Objective_function();
  MCMC_output();
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::Convert_log_parameters(void)
{
  ofstream& report1= *pad_report1;
 // Assume close to equilibrium at f=0 at the start   
 //Use the initial recruitment dev to fill out the initial age composition
  initN(rcrage+1) = mfexp(mean_log_recruit +  dev_log_recruit(styr) - M(rcrage));
  for (j=rcrage+2;j<=trmage;j++) {
    initN(j) = initN(j-1)*mfexp(-M(j));
  }
  initN(trmage) /= (1.0 - mfexp(-M(trmage)));
 //devs for initN are turned off 
 for (j=rcrage+1;j<=trmage;j++) {
   initN(j) = initN(j)*mfexp(dev_log_initN(j));
 }
 log_recruit=mean_log_recruit +  dev_log_recruit;
 recruit = mfexp(log_recruit);
 // Acoustic survey random walk setup
 for (i=styr;i<=endyr;i++) {
   log_q1(i)=log_q1_mean+log_q1_dev(i);
   log_q2(i)=log_q2_mean+log_q2_dev(i);
   log_q3(i)=log_q3_mean+log_q3_dev(i); 
   q1(i)=mfexp(log_q1(i));
   q2(i)=mfexp(log_q2(i));
   q3(i)=mfexp(log_q3(i));
 }
 F = mfexp(mean_log_F + dev_log_F);
 q4 = mfexp(log_q4);
 q5 = mfexp(log_q5);
 q6 = mfexp(log_q6);
}

void model_parameters::Selectivity(void)
{
  ofstream& report1= *pad_report1;
 // Fishery selectivity
 for (i=styr;i<=endyr;i++) {
   slp1_fsh(i)=mfexp(log_slp1_fsh_mean+slp1_fsh_dev(i));
   inf1_fsh(i)=inf1_fsh_mean+inf1_fsh_dev(i);
   slp2_fsh(i)=mfexp(log_slp2_fsh_mean+slp2_fsh_dev(i));
   inf2_fsh(i)=inf2_fsh_mean+inf2_fsh_dev(i);
   for (j=rcrage;j<=trmage;j++) {
     slctfsh(i,j) = (1/(1+mfexp(-(slp1_fsh(i))*(double(j)-(inf1_fsh(i))))))*
       (1-1/(1+mfexp(-(slp2_fsh(i))*(double(j)-(inf2_fsh(i))))));
   }
 // The plan would be to check and adjust the max selected age as needed
 //if (i<1974)
   //slctfsh(i)=slctfsh(i)/slctfsh(i,3);
 //else
   //slctfsh(i)=slctfsh(i)/slctfsh(i,7);
   //slctfsh(i)=slctfsh(i)/mean(slctfsh(i)(4,8));
 }
 M(1)=0.9; M(2)=0.45; M(3)=0.3; M(4)=0.3; M(5)=0.3;
 M(6)=0.30; M(7)=0.30; M(8)=0.30; M(9)=0.30; M(10)=0.30;
 //0.9	0.45	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3	0.3
 for(int i=1; i<=10; i++) M(i)*=natMscalar;
 for (j=rcrage;j<=trmage;j++) {
   slctsrv1(j) = (1-1/(1+mfexp(-mfexp(log_slp2_srv1)*(double(j)-inf2_srv1))));
 }
 slctsrv1=slctsrv1/slctsrv1(3);
 slctsrv1(rcrage)=0;
 slctsrv1(rcrage+1)=0;	
 for (j=rcrage;j<=trmage;j++) {
   slctsrv2(j) = (1/(1+mfexp(-mfexp(log_slp1_srv2)*(double(j)-inf1_srv2))))
     *(1-1/(1+mfexp(-mfexp(log_slp2_srv2)*(double(j)-inf2_srv2))));
 }
 slctsrv2=slctsrv2/slctsrv2(10);
  for (j=rcrage;j<=trmage;j++) {
    slctsrv3(j) = (1/(1+mfexp(-mfexp(log_slp1_srv3)*(double(j)-inf1_srv3))));
  }
  slctsrv3=slctsrv3/slctsrv3(10);
 // Survey 6 selectivity
  for (j=rcrage;j<=trmage;j++) {
    slctsrv6(j) = (1/(1+mfexp(-mfexp(log_slp1_srv6)*(double(j)-inf1_srv6))))
      *(1-1/(1+mfexp(-mfexp(log_slp2_srv6)*(double(j)-inf2_srv6))));
  }
  slctsrv6=slctsrv6/slctsrv6(1);
  // Calculate uncertainty in logit space where it makes more
  // sense. But it does break when selex=0 or 1 identically so
  // need to be careful there.
  slctsrv1_logit = -log(1/(slctsrv1-1e-10)-1);
  slctsrv2_logit = -log(1/(slctsrv2-1e-10)-1);
  slctsrv3_logit = -log(1/(slctsrv3-1e-10)-1);
  slctsrv6_logit = -log(1/(slctsrv6-1e-10)-1);
  slctfsh_logit  = -log(1/(slctfsh(endyr-1)-1e-10)-1);
}

void model_parameters::Mortality(void)
{
  ofstream& report1= *pad_report1;
  for (i=styr;i<=endyr;i++) {
    for (j=rcrage;j<=trmage;j++) {
      Z(i,j)=(F(i)*slctfsh(i,j))+M(j);
    }
  }  
}

void model_parameters::Numbers_at_age(void)
{
  ofstream& report1= *pad_report1;
  N(styr)(rcrage+1,trmage)=initN;
  for (i=styr;i<=endyr;i++) {
    N(i,rcrage)=recruit(i);
  }
  for (i=styr;i<endyr;i++) {
    for (j=rcrage;j<trmage;j++) {
      N(i+1,j+1)=N(i,j)*mfexp(-Z(i,j));
    }  
    N(i+1,trmage)+=N(i,trmage)*mfexp(-Z(i,trmage));
  }
  endN=N(endyr);
}

void model_parameters::Catch_at_age(void)
{
  ofstream& report1= *pad_report1;
 // Catch at age and survey numbers at age
  for (i=styr;i<=endyr;i++) {
    for (j=rcrage;j<=trmage;j++) {
      C(i,j)=N(i,j)*((F(i)*slctfsh(i,j))/Z(i,j))*(1-mfexp(-Z(i,j)));
      Eec(i,j)=N(i,j)*(M(j)/Z(i,j))*(1-mfexp(-Z(i,j)));
      Nsrv1(i,j)=slctsrv1(j)*N(i,j)*mfexp(-yrfrct_srv1(i)*Z(i,j));
      Nsrv2(i,j)=slctsrv2(j)*N(i,j)*mfexp(-yrfrct_srv2(i)*Z(i,j));
      Nsrv3(i,j)=slctsrv3(j)*N(i,j)*mfexp(-yrfrct_srv3(i)*Z(i,j));
      Nsrv6(i,j)=slctsrv6(j)*N(i,j)*mfexp(-yrfrct_srv6(i)*Z(i,j));
     }
  }
}

void model_parameters::Expected_values(void)
{
  ofstream& report1= *pad_report1;
 for (i=styr;i<=endyr;i++){
   Ecattot(i) = 1000000*sum(elem_prod(C(i),wt_fsh(i)));
   Eecocon(i) = 1000000*sum(elem_prod(Eec(i),wt_pop(i)));
   Ecatp(i) = (C(i)/sum(C(i)))*age_trans;
   Elenp(i) = Ecatp(i) * len_trans1;
   Eindxsurv1(i)= q1(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv1(i)*Z(i))),slctsrv1),wt_srv1(i)));
   Esrvp1(i) = (Nsrv1(i)/sum(Nsrv1(i)))*age_trans;
   Esrvlenp1(i) = Esrvp1(i) * len_trans3;
   Eindxsurv2(i)= q2(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv2(i)*Z(i))),slctsrv2),wt_srv2(i)));
   Esrvp2(i) = (Nsrv2(i)/sum(Nsrv2(i)))*age_trans;
   Esrvlenp2(i) = Esrvp2(i) * len_trans2;
   Eindxsurv3(i)= q3(i)*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv3(i)*Z(i))),slctsrv3),wt_srv3(i)));
   Esrvp3(i) = (Nsrv3(i)/sum(Nsrv3(i)))*age_trans;
   Esrvlenp3(i) = Esrvp3(i) * len_trans2;
   Eindxsurv4(i)= q4*pow(N(i,1),(q4_pow+1));
   Eindxsurv5(i)= q5*pow(N(i,2),(q5_pow+1));	
   Eindxsurv6(i)= q6*sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-yrfrct_srv6(i)*Z(i))),slctsrv6),wt_srv6(i)));
   Esrvp6(i) = (Nsrv6(i)/sum(Nsrv6(i)))*age_trans;
   Esrvlenp6(i) = Esrvp6(i) * len_trans2;
   // 3+ biomass
   Esumbio(i)= N(i)(rcrage+2,trmage)*wt_pop(i)(rcrage+2,trmage);
   // Total biomass
   Etotalbio(i)= N(i)(rcrage,trmage)*wt_pop(i)(rcrage,trmage);
   // 2+ biomass at spawning (use for apportioning to management area)
   Espawnbio_2plus(i)= sum(elem_prod(elem_prod(N(i)(rcrage+1,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage+1,trmage))),wt_srv1(i)(rcrage+1,trmage)));
   // //1+ biomass at spawning
   //     Esumbio(i)= sum(elem_prod(elem_prod(N(i)(rcrage,trmage),mfexp(-yrfrct_srv1(i)*Z(i)(rcrage,trmage))),wt_srv1(i)(rcrage,trmage)));
   Espawnbio(i)= sum(elem_prod(elem_prod(elem_prod(N(i),mfexp(-0.21*Z(i))),wt_spawn(i)),0.5*mat));
 }
 Espawnbio_log=log(Espawnbio);
 Esumbio_log=log(Esumbio);
 // Do upper and lower accumulation in expected age composition
 // Fishery
  for (i=1;i<=nyrs_fsh;i++) {
    if(fshyrs(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_fsh(i)) {
	Ecatp(fshyrs(i),ac_yng_fsh(i)) += Ecatp(fshyrs(i),j);
	Ecatp(fshyrs(i),j) = 0;
      }
      if(j>ac_old_fsh(i)) {
	Ecatp(fshyrs(i),ac_old_fsh(i)) += Ecatp(fshyrs(i),j);
	Ecatp(fshyrs(i),j) = 0;
      }
    }
  } 
 // Survey 1
  for (i=1;i<=nyrsac_srv1;i++) {
    if(srv_acyrs1(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv1(i)) {
	Esrvp1(srv_acyrs1(i),ac_yng_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
	Esrvp1(srv_acyrs1(i),j) = 0;
      }
      if(j>ac_old_srv1(i)) {
	Esrvp1(srv_acyrs1(i),ac_old_srv1(i)) += Esrvp1(srv_acyrs1(i),j);
	Esrvp1(srv_acyrs1(i),j) = 0;
      }
    }
  }
  // Survey 2
  for (i=1;i<=nyrsac_srv2;i++) {
    if(srv_acyrs2(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv2(i)) {
	Esrvp2(srv_acyrs2(i),ac_yng_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
	Esrvp2(srv_acyrs2(i),j) = 0.;
      }
      if(j>ac_old_srv2(i)) {
	Esrvp2(srv_acyrs2(i),ac_old_srv2(i)) += Esrvp2(srv_acyrs2(i),j);
	Esrvp2(srv_acyrs2(i),j) = 0;
      }
    }
  }
  // Survey 6
  for (i=1;i<=nyrsac_srv6;i++) {
    if(srv_acyrs6(i)>endyr) break;
    for (j=rcrage;j<=trmage;j++) {
      if(j<ac_yng_srv6(i)) {
	Esrvp2(srv_acyrs6(i),ac_yng_srv6(i)) += Esrvp6(srv_acyrs6(i),j);
	Esrvp6(srv_acyrs6(i),j) = 0.;
      }
      if(j>ac_old_srv6(i)) {
	Esrvp6(srv_acyrs6(i),ac_old_srv6(i)) += Esrvp6(srv_acyrs6(i),j);
	Esrvp6(srv_acyrs6(i),j) = 0;
      }
    }
  }
}

void model_parameters::Projections(void)
{
  ofstream& report1= *pad_report1;
 for (i=endyr+1;i<=endyr+5;i++) {
   //  recruitment with bias correction
   //    recruit_proj(i)=mfexp(log_recr_proj(i)+(sigmasq_recr/2));
   // for MCMC projections to get the prob < B20
   //   recruit_proj(i)=mfexp(log_recr_proj(i));
   // or just use average recruitment after 1977
   // note that endyr-1 is the last year for mean
   // standard method
   recruit_proj(i)=mean(recruit(1978,endyr-1));
 }
 for (i=endyr+1;i<=endyr+5;i++) {
   N_proj(i,rcrage)=recruit_proj(i);
 }
 for (j=rcrage;j<trmage;j++) {
   N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
 }  
 N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));
 // Set 2007 year class to mean
 // note that endyr-1 is the last year for mean
 //   for (j=rcrage;j<trmage;j++)
 //   {
 //   N_proj(endyr+1,j+1)=N(endyr,j)*mfexp(-Z(endyr,j));
 //   }  
 //   N_proj(endyr+1,rcrage+1)=mean(recruit(1979,endyr-1))*mfexp(-Z(endyr,rcrage));
 //   N_proj(endyr+1,trmage)+=N(endyr,trmage)*mfexp(-Z(endyr,trmage));
 // Averaging window for selectivity 
 endyr_avg_slct=endyr-1;
 styr_avg_slct=endyr_avg_slct-4;
 for (j=rcrage;j<=trmage;j++) {
   slctfsh_proj(j) = 0;
   for (i=styr_avg_slct;i<=endyr_avg_slct;i++) {
     slctfsh_proj(j) += slctfsh(i,j);
   }
 }
 slctfsh_proj=slctfsh_proj/max(slctfsh_proj);
 //Forward projections
  for (i=endyr+1;i<=endyr+5;i++) {
    // have to get Z twice
    //  Tuning loop to get the spawning biomass adjustment right
    F_proj(i)=Ftarget(i+retro_yrs);
    for (loop=1;loop<=20;loop++) {
      for (j=rcrage;j<=trmage;j++) {
	Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M(j);
      }  
      sbio = sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
      //  Set the fishing mortality rate
      F_proj(i)=Ftarget(i+retro_yrs);
      if (sbio < B40) {
	F_proj(i)=Ftarget(i+retro_yrs)*(((sbio/B40)-0.05)/(1-0.05));
	// SSL control rule
	//   F_proj(i)=Ftarget(i+retro_yrs)*(((sbio/B40)-0.2)/(1-0.2));
      }
    }
  // Total mortality
    for (j=rcrage;j<=trmage;j++) {
      Z_proj(i,j)=(F_proj(i)*slctfsh_proj(j))+M(j);
    }  
  //  Numbers at age
   if(i<endyr+5) {
     for (j=rcrage;j<trmage;j++) {
       N_proj(i+1,j+1)=N_proj(i,j)*mfexp(-Z_proj(i,j));
     }  
     N_proj(i+1,trmage)+=N_proj(i,trmage)*mfexp(-Z_proj(i,trmage)); 
   }
  // Catches
  for (j=rcrage;j<=trmage;j++) {
    C_proj(i,j)=N_proj(i,j)*((F_proj(i)*slctfsh_proj(j))/Z_proj(i,j))*(1-mfexp(-Z_proj(i,j)));
    //    Nsrv_proj(i,j)=q2*slctsrv1(j)*N_proj(i,j)*mfexp(-yrfrct_srv1(endyr)*Z_proj(i,j));  
    Nsrv_proj(i,j)=N_proj(i,j)*mfexp(-yrfrct_srv6(endyr)*Z_proj(i,j));  
  }
  //  Total catches and biomass
  Ecattot_proj(i) = 1000000*sum(elem_prod(C_proj(i),wt_fsh_proj));
  // 3+ biomass
  Esumbio_proj(i)= N_proj(i)(rcrage+2,trmage)*wt_pop_proj(rcrage+2,trmage);
  // Alternative: 2+ biomass
  //    Esumbio_proj(i)= N_proj(i)(rcrage+1,trmage)*wt_pop_proj(rcrage+1,trmage);
  Exrate_proj(i)=Ecattot_proj(i)/(1000000*Esumbio_proj(i));
  Espawnbio_proj(i)= sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-0.21*Z_proj(i))),wt_spawn_proj),0.5*mat));
  //Summer acoustic
  //    Esrv_proj(i)= q6*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv6(endyr)*Z_proj(i))),slctsrv6),wt_srv_proj));
  //Winter acoutic 
  Esrv_proj(i)= q1(endyr)*sum(elem_prod(elem_prod(elem_prod(N_proj(i),mfexp(-yrfrct_srv1(endyr)*Z_proj(i))),slctsrv1),wt_srv_proj));
  }
}

void model_parameters::Objective_function(void)
{
  ofstream& report1= *pad_report1;
 // Fishery likelihoods
 //Total catch
  //loglik(1) = -.5*norm2(elem_div((log(cattot)-log(Ecattot)),cattot_log_sd));
 loglik(1)=0;
 for(i=styr; i<=endyr;i++){
   if(i>endyr) break;
   loglik(1) += -.5*square((log(cattot(i))-log(Ecattot(i)))/cattot_log_sd(i));
 }	  
 loglik(2)=0;
 for (i=1;i<=nyrs_fsh;i++) {
    if(fshyrs(i)>endyr) break; 	// ignore data after retroyear
    llcatp(i) = 0;
    for (j=ac_yng_fsh(i);j<=ac_old_fsh(i);j++) {
      llcatp(i) += multN_fsh(i)*(catp(i,j)+o)*log((Ecatp(fshyrs(i),j)+o)/(catp(i,j)+o));
      res_fish(i,j)=catp(i,j);
      res_fish(i,trmage-rcrage+j+1)=Ecatp(fshyrs(i),j);
      if(multN_fsh(i)>0) {
	// pearson_fish(i,j)=(catp(i,j)-Ecatp(fshyrs(i),j))/sqrt((Ecatp(fshyrs(i),j)*(1.-Ecatp(fshyrs(i),j)))/multN_fsh(i));	
      }
    }
    if(multN_fsh(i)>0) {	
      effN_fsh(i) = sum(elem_prod(Ecatp(fshyrs(i)),(1-Ecatp(fshyrs(i)))))/sum(square(catp(i)-Ecatp(fshyrs(i))));
    }
    loglik(2) += llcatp(i);
 }
 //Length composition
 loglik(3)=0;
 for (i=1;i<=nyrslen_fsh;i++) {
   if(fshlenyrs(i)>endyr) break;
   lllenp(i) = 0;
   for (j=1;j<=nbins1;j++) {
     lllenp(i) += multNlen_fsh(i)*(lenp(i,j)+o)*log((Elenp(fshlenyrs(i),j)+o)/(lenp(i,j)+o));
   }
   loglik(3) += lllenp(i);
 }
 // Survey 1 likelihoods
 // Total biomass  
  loglik(4)=0;
  for(i=1; i<=nyrs_srv1;i++){
    if(srvyrs1(i)>endyr) break;
    loglik(4)+=-.5*square((log(indxsurv1(i))-log(Eindxsurv1(srvyrs1(i)))+square(indxsurv_log_sd1(i))/2.)/indxsurv_log_sd1(i));
  }
  RMSE_srv1=0;
  if(!isretro)
    RMSE_srv1= sqrt(norm2(log(indxsurv1)-log(Eindxsurv1(srvyrs1))+square(indxsurv_log_sd1)/2.)/nyrs_srv1);
 //Age composition
  loglik(5)=0;
  for (i=1;i<=nyrsac_srv1;i++) {
    if(srv_acyrs1(i)>endyr) break;
    llsrvp1(i) = 0;
    for (j=ac_yng_srv1(i);j<=ac_old_srv1(i);j++) {
      llsrvp1(i) += multN_srv1(i)*(srvp1(i,j)+o)*log((Esrvp1(srv_acyrs1(i),j)+o)/(srvp1(i,j)+o));
      res_srv1(i,j)=srvp1(i,j);
      res_srv1(i,trmage-rcrage+j+1)=Esrvp1(srv_acyrs1(i),j);
      if(multN_srv1(i)>0) {
	// pearson_srv1(i,j)=(srvp1(i,j)-Esrvp1(srv_acyrs1(i),j))/sqrt((Esrvp1(srv_acyrs1(i),j)*(1.-Esrvp1(srv_acyrs1(i),j)))/multN_srv1(i));	
      }
    }
    if(multN_srv1(i)>0) {
      effN_srv1(i) = sum(elem_prod(Esrvp1(srv_acyrs1(i)),(1-Esrvp1(srv_acyrs1(i)))))/sum(square(srvp1(i)-Esrvp1(srv_acyrs1(i))));
    }
    loglik(5) += llsrvp1(i);
  }
 //length composition
 loglik(6)=0;
 for (i=1;i<=nyrslen_srv1;i++) {
   if(srv_lenyrs1(i)>endyr) break;
   llsrvlenp1(i) = 0;
   for (j=1;j<=nbins3;j++) {
     llsrvlenp1(i) += multNlen_srv1(i)*(srvlenp1(i,j)+o)*log((Esrvlenp1(srv_lenyrs1(i),j)+o)/(srvlenp1(i,j)+o));
   }
   loglik(6) += llsrvlenp1(i);
 }
 // Survey 2 likelihoods
 //Total biomass    
  loglik(7) =0;
  for (i=1;i<=nyrs_srv2;i++){
    if(srvyrs2(i)>endyr) break;
    loglik(7)+=-.5*square((log(indxsurv2(i))-log(Eindxsurv2(srvyrs2(i)))+square(indxsurv_log_sd2(i))/2.)/indxsurv_log_sd2(i));
  }
  RMSE_srv2=0;
  if(!isretro)
    RMSE_srv2= sqrt(norm2(log(indxsurv2)-log(Eindxsurv2(srvyrs2))+square(indxsurv_log_sd2)/2.)/nyrs_srv2);
 // Age composition
  loglik(8)=0;
  for (i=1;i<=nyrsac_srv2;i++) {
    llsrvp2(i) = 0;
   if(srv_acyrs2(i)>endyr) break;
   for (j=ac_yng_srv2(i);j<=ac_old_srv2(i);j++) {
     llsrvp2(i) += multN_srv2(i)*(srvp2(i,j)+o)*log((Esrvp2(srv_acyrs2(i),j)+o)/(srvp2(i,j)+o));
     res_srv2(i,j)=srvp2(i,j);
     res_srv2(i,trmage-rcrage+j+1)=Esrvp2(srv_acyrs2(i),j);
     if(multN_srv2(i)>0) {
       // pearson_srv2(i,j)=(srvp2(i,j)-Esrvp2(srv_acyrs2(i),j))/sqrt((Esrvp2(srv_acyrs2(i),j)*(1.-Esrvp2(srv_acyrs2(i),j)))/multN_srv2(i));	
     }  
   }
   if(multN_srv2(i)>0) {	
     effN_srv2(i) = sum(elem_prod(Esrvp2(srv_acyrs2(i)),(1-Esrvp2(srv_acyrs2(i)))))/sum(square(srvp2(i)-Esrvp2(srv_acyrs2(i))));
   }
   loglik(8) += llsrvp2(i);
  }
  // length composition
  loglik(9)=0;
  for (i=1;i<=nyrslen_srv2;i++) {
    if(srv_lenyrs2(i) > endyr) break;
    llsrvlenp2(i) = 0;
    for (j=1;j<=nbins2;j++) {
      llsrvlenp2(i) += multNlen_srv2(i)*(srvlenp2(i,j)+o)*log((Esrvlenp2(srv_lenyrs2(i),j)+o)/(srvlenp2(i,j)+o));
    }
    loglik(9) += llsrvlenp2(i);
  }
  loglik(10) = 0;
    loglik(11)=0;
    for(i=1; i<=nyrs_srv3;i++){
      if(srvyrs3(i)>endyr) break;
      loglik(11) += -.5*square((log(indxsurv3(i))-log(Eindxsurv3(srvyrs3(i)))+square(indxsurv_log_sd3(i))/2.)/indxsurv_log_sd3(i));
    }
    RMSE_srv3=0;
    if(!isretro)
      RMSE_srv3= sqrt(norm2(log(indxsurv3)-log(Eindxsurv3(srvyrs3))+square(indxsurv_log_sd3)/2.)/nyrs_srv3);
 // age composition
    loglik(12)=0;
  for (i=1;i<=nyrsac_srv3;i++) {
    if(srv_acyrs3(i)>endyr) break;
    llsrvp3(i) = 0;
    for (j=rcrage;j<=trmage;j++) {
      llsrvp3(i) += multN_srv3(i)*(srvp3(i,j)+o)*log((Esrvp3(srv_acyrs3(i),j)+o)/(srvp3(i,j)+o));
      res_srv3(i,j)=srvp3(i,j);
      res_srv3(i,trmage-rcrage+j+1)=Esrvp3(srv_acyrs3(i),j);
      if(multN_srv3(i)>0) {
	// pearson_srv3(i,j)=(srvp3(i,j)-Esrvp3(srv_acyrs3(i),j))/sqrt((Esrvp3(srv_acyrs3(i),j)*(1.-Esrvp3(srv_acyrs3(i),j)))/multN_srv3(i));	
      }
    }
    if(multN_srv3(i)>0) {		
      effN_srv3(i) = sum(elem_prod(Esrvp3(srv_acyrs3(i)),(1-Esrvp3(srv_acyrs3(i)))))/sum(square(srvp3(i)-Esrvp3(srv_acyrs3(i))));	
    }
    loglik(12) += llsrvp3(i);
  }
 // length composition
  loglik(13)=0;
  for (i=1;i<=nyrslen_srv3;i++) {
    if(srv_lenyrs3(i)>endyr) break;
    llsrvlenp3(i) = 0;
    for (j=1;j<=nbins2;j++) {
      llsrvlenp3(i) += multNlen_srv3(i)*(srvlenp3(i,j)+o)*log((Esrvlenp3(srv_lenyrs3(i),j)+o)/(srvlenp3(i,j)+o));
      res_srv3len(i,j)=srvlenp3(i,j);
      res_srv3len(i,nbins2+j)=Esrvlenp3(srv_lenyrs3(i),j);
      if(multNlen_srv3(i)>0) {
	// pearson_srv3len(i,j)=(srvlenp3(i,j)-Esrvlenp3(srv_lenyrs3(i),j))/sqrt((Esrvlenp3(srv_lenyrs3(i),j)*(1.-Esrvlenp3(srv_lenyrs3(i),j)))/multNlen_srv3(i));	
      }
    }
    loglik(13) += llsrvlenp3(i);
  }
  loglik(14)=0; loglik(15)=0;
  for(i=1; i<=nyrs_srv4;i++){ 	// assuming srv4 and srv5 have identical structure
    if(srvyrs4(i) >endyr) break;
    loglik(14) += -.5*square((log(indxsurv4(i))-log(Eindxsurv4(srvyrs4(i)))+square(indxsurv_log_sd4(i))/2.)/indxsurv_log_sd4(i));
   //   loglik(14) = 0;
    loglik(15) += -.5*square((log(indxsurv5(i))-log(Eindxsurv5(srvyrs5(i)))+square(indxsurv_log_sd5(i))/2.)/indxsurv_log_sd5(i));
  }
   RMSE_srv5=0;
   if(!isretro)
     RMSE_srv5= sqrt(norm2(log(indxsurv5)-log(Eindxsurv5(srvyrs5))+square(indxsurv_log_sd5)/2.)/nyrs_srv5);
   RMSE_srv4=0;
   if(!isretro)
     RMSE_srv4= sqrt(norm2(log(indxsurv4)-log(Eindxsurv4(srvyrs4))+square(indxsurv_log_sd4)/2.)/nyrs_srv4);
  //Total biomass    
 loglik(16)=0;
 for(i=1;i<=nyrs_srv6;i++){
   if(srvyrs6(i)>endyr) break;
   loglik(16)+=-.5*square((log(indxsurv6(i))-log(Eindxsurv6(srvyrs6(i)))+square(indxsurv_log_sd6(i))/2.)/indxsurv_log_sd6(i));
 }
 RMSE_srv6=0;
 if(!isretro)
   RMSE_srv6= sqrt(norm2(log(indxsurv6)-log(Eindxsurv6(srvyrs6))+square(indxsurv_log_sd6)/2.)/nyrs_srv6);
 // Age composition
  loglik(17)=0;
 for (i=1;i<=nyrsac_srv6;i++) {
   if(srv_acyrs6(i)>endyr) break;
   llsrvp6(i) = 0;
   for (j=ac_yng_srv6(i);j<=ac_old_srv6(i);j++) {
     llsrvp6(i) += multN_srv6(i)*(srvp6(i,j)+o)*log((Esrvp6(srv_acyrs6(i),j)+o)/(srvp6(i,j)+o));
     res_srv6(i,j)=srvp6(i,j);
     res_srv6(i,trmage-rcrage+j+1)=Esrvp6(srv_acyrs6(i),j);
     if(multN_srv6(i)>0) {
       // pearson_srv6(i,j)=(srvp6(i,j)-Esrvp6(srv_acyrs6(i),j))/sqrt((Esrvp6(srv_acyrs6(i),j)*(1.-Esrvp6(srv_acyrs6(i),j)))/multN_srv6(i));	
     }  
   }
   if(multN_srv6(i)>0) {	
     effN_srv6(i) = sum(elem_prod(Esrvp6(srv_acyrs6(i)),(1-Esrvp6(srv_acyrs6(i)))))/sum(square(srvp6(i)-Esrvp6(srv_acyrs6(i))));
   }
   loglik(17) +=llsrvp6(i);
 }
  for (i=1;i<=nyrslen_srv6;i++) {
    if(srv_lenyrs6(i)>endyr) break;
    llsrvlenp6(i) = 0;
  for (j=1;j<=nbins2;j++) {
    llsrvlenp6(i) += multNlen_srv6(i)*(srvlenp6(i,j)+o)*log((Esrvlenp6(srv_lenyrs6(i),j)+o)/(srvlenp6(i,j)+o));
  }
  loglik(17) += llsrvlenp6(i);
  }
  // Constraints on recruitment. Assumed sigmaR=1.3 for all devs
  loglik(18)= 0;
  loglik(18) += -0.5*norm2(dev_log_recruit/sigmaR);
 // Normal process error on selectivity deviations. Note
 // rwlk_sd(styr,endyr-1) b/c if using retro they will be too
 // long since read in as data with the original endyr value
  loglik(19)  = -0.5*norm2(elem_div(first_difference(slp1_fsh_dev),rwlk_sd(styr,endyr-1))); 
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf1_fsh_dev),4.0*rwlk_sd(styr,endyr-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(slp2_fsh_dev),rwlk_sd(styr,endyr-1)));
  loglik(19) += -0.5*norm2(elem_div(first_difference(inf2_fsh_dev),4.0*rwlk_sd(styr,endyr-1)));
 // Recruitment in projection mode, but skip it if doing retro
 // peels b/c it makes no sense and breaks this code
  if(last_phase()) {
    loglik(20) =  -(1/(2.0*sigmasq_recr))*norm2(log_recr_proj - log_mean_recr_proj);
  } else {
    loglik(20)=0;
  }
 // Normal process error on catchability deviations. Note
 // rwlk_sd(styr,endyr-1) b/c if using retro they will be too
 // long since read in as data with the original endyr value
  loglik(21)  = -0.5*norm2(elem_div(first_difference(log_q1_dev),q1_rwlk_sd(styr,endyr-1))); 
  loglik(21)  += -0.5*norm2(elem_div(first_difference(log_q2_dev),q2_rwlk_sd(styr,endyr-1))); 
  loglik(21)  += -0.5*norm2(elem_div(first_difference(log_q3_dev),q3_rwlk_sd(styr,endyr-1))); 
  loglik(22)= 0;
 // Prior on trawl catchability       
 loglik(23) = -.5*square((log_q2_mean-log(0.85))/0.1);
 loglik(24) = 0;
 // these broad priors stabilize estimation, particularly for
 // retros, and imply uniform 1 selex for this survey
 loglik(24) -= dnorm(log_slp2_srv6, 0,2, true);
 loglik(24) -= dnorm(inf2_srv6, 10,3, true);
 objfun = -sum(loglik);
 // Variable to do a likelihood profile over
  var_prof=Esumbio(endyr);
}

void model_parameters::MCMC_output(void)
{
  ofstream& report1= *pad_report1;
 if(mceval_phase()) {
   // Dumping data to mceval.dat
   report1<<mean(recruit(1978,endyr-1));
   report1<<" ";
   report1<<endyr-2;
   report1<<" ";
   report1<<F(endyr-2);
   report1<<" ";
   report1<<Espawnbio(endyr-2);
   report1<<" ";
   report1<<recruit(endyr-2);
   report1<<" ";
   report1<<endyr-1;
   report1<<" ";
   report1<<F(endyr-1);
   report1<<" ";
   report1<<Espawnbio(endyr-1);
   report1<<" ";
   report1<<recruit(endyr-1);
   report1<<" ";
   report1<<endyr;
   report1<<" ";
   report1<<F(endyr);
   report1<<" ";
   report1<<Espawnbio(endyr);
   report1<<" ";
   report1<<recruit(endyr);
   report1<<" ";
   report1<<endyr+1;
   report1<<" ";
   report1<<F_proj(endyr+1);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+1);
   report1<<" ";
   report1<<recruit_proj(endyr+1);
   report1<<" ";
   report1<<endyr+2;
   report1<<" ";
   report1<<F_proj(endyr+2);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+2);
   report1<<" ";
   report1<<recruit_proj(endyr+2);
   report1<<" ";
   report1<<endyr+3;
   report1<<" ";
   report1<<F_proj(endyr+3);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+3);
   report1<<" ";
   report1<<recruit_proj(endyr+3);
   report1<<" ";
   report1<<endyr+4;
   report1<<" ";
   report1<<F_proj(endyr+4);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+4);
   report1<<" ";
   report1<<recruit_proj(endyr+4);
   report1<<" ";
   report1<<endyr+5;
   report1<<" ";
   report1<<F_proj(endyr+5);
   report1<<" ";
   report1<<Espawnbio_proj(endyr+5);
   report1<<" ";
   report1<<recruit_proj(endyr+5);
   report1 << endl;
  }
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e0, 1.e-1, 1.e-4, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100, 1000, 1000, 100000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Objective function" << endl << objfun << endl;
  report << "Likelihood components" << endl << loglik << endl;
  report << "Recruits" << endl <<  recruit << endl; 
  report << "Expected spawning biomass" << endl << Espawnbio << endl;
  report << "Expected total biomass" << endl << Etotalbio << endl;
  report << "Expected summary (age 3+) biomass" << endl << Esumbio << endl;
  report << "Expected spawning biomass age 2+" << endl << Espawnbio_2plus << endl;
  report << "Fishing mortalities" << endl <<  F << endl ; 
  report << "Ecosystem comsumption" << endl <<  Eecocon << endl; 
  report << "Initial age comp" << endl <<  initN << endl; 
  report << "Total catch" << endl <<  cattot << endl;  
  report << "Expected total catch" << endl <<  Ecattot << endl;  
  report << "Natural mortality" << endl <<  M << endl; 
  report << "Numbers at age" << endl << N << endl << endl;
  report << "Fishery selectivity means" << endl;  
  report << log_slp1_fsh_mean << " " <<  inf1_fsh_mean << " " << log_slp2_fsh_mean  << " " << inf2_fsh_mean << endl; 
  report << "Fishery selectivity deviances" << endl;  
  report << slp1_fsh_dev << endl; 
  report << inf1_fsh_dev << endl; 
  report << slp2_fsh_dev << endl; 
  report << inf2_fsh_dev << endl; 
  report << "Fishery selectivity vectors" << endl;  
  report <<  slp1_fsh << endl; 
  report <<  inf1_fsh << endl; 
  report <<  slp2_fsh << endl; 
  report <<  inf2_fsh << endl; 
  report << "Fishery selectivity" << endl << slctfsh << endl;
  report << "Fishery age composition likelihoods" << endl << llcatp << endl;
  report << "Fishery age composition" << endl << catp << endl;
  report << "Fishery expected age composition" << endl << Ecatp << endl;
  report << "Fishery observed and expected age comp" << endl << res_fish << endl;
  report << "Fishery Pearson residuals age comp" << endl << pearson_fish << endl; 
  report << "Fishery input N" << endl << multN_fsh << endl;  
  report << "Fishery effective N age comp" << endl << effN_fsh << endl;   
  report << "Fishery length composition likelihoods" << endl << lllenp << endl;
  report << "Fishery length composition" << endl << lenp << endl;
  report << "Fishery expected length composition" << endl << Elenp << endl << endl;
  report << "Survey 1 q" << endl << q1 << endl;
  report << "Survey 1 selectivity parameters" << endl;  
  report << log_slp2_srv1 << " " << inf2_srv1 << endl; 
  report << "Survey 1 selectivity" << endl << slctsrv1 << endl;
  report << "Survey 1 expected index" << endl << Eindxsurv1 << endl;
  report << "Survey 1 RMSE" << endl << RMSE_srv1 << endl;
  report << "Survey 1 age composition likelihoods" << endl << llsrvp1 << endl;
  report << "Survey 1 age composition" << endl << srvp1 << endl;
  report << "Survey 1 expected age composition" << endl << Esrvp1 << endl;
  report << "Survey 1 observed and expected age comp" << endl << res_srv1 << endl;
  report << "Survey 1 Pearson residuals age comp" << endl << pearson_srv1 << endl; 
  report << "Survey 1 input N" << endl << multN_srv1 << endl;  
  report << "Survey 1 effective N age comp" << endl << effN_srv1 << endl;   
  report << "Survey 1 length composition likelihoods" << endl << llsrvlenp1 << endl;
  report << "Survey 1 length composition" << endl << srvlenp1 << endl;
  report << "Survey 1 expected length composition" << endl << Esrvlenp1 << endl << endl;
  report << "Survey 2 q" << endl << q2 << endl;
  report << "Survey 2 selectivity parameters" << endl;  
  report << log_slp1_srv2 << " " << inf1_srv2 << " " << log_slp2_srv2 << " " << inf2_srv2 << endl; 
  report << "Survey 2 selectivity" << endl << slctsrv2 << endl;
  report << "Survey 2 expected index" << endl << Eindxsurv2 << endl;
  report << "Survey 2 RMSE" << endl << RMSE_srv2 << endl;
  report << "Survey 2 age composition likelihoods" << endl << llsrvp2 << endl;
  report << "Survey 2 age composition" << endl << srvp2 << endl;
  report << "Survey 2 expected age composition" << endl << Esrvp2 << endl;
  report << "Survey 2 observed and expected age comp" << endl << res_srv2 << endl;
  report << "Survey 2 Pearson residuals age comp" << endl << pearson_srv2 << endl;  
  report << "Survey 2 input N" << endl << multN_srv2 << endl;  
  report << "Survey 2 effective N age comp" << endl << effN_srv2 << endl;   
  report << "Survey 2 length composition likelihoods" << endl << llsrvlenp2 << endl;
  report << "Survey 2 length composition" << endl << srvlenp2 << endl;
  report << "Survey 2 expected length composition" << endl << Esrvlenp2 << endl << endl;
  report << "Survey 3 q" << endl << q3 << endl;
  report << "Survey 3 selectivity parameters" << endl;  
  report << log_slp1_srv3 << " " <<  inf1_srv3 << endl; 
  report << "Survey 3 selectivity" << endl << slctsrv3 << endl;
  report << "Survey 3 expected index" << endl << Eindxsurv3 << endl;
  report << "Survey 3 RMSE" << endl << RMSE_srv3 << endl;
  report << "Survey 3 age composition likelihoods" << endl << llsrvp3 << endl;
  report << "Survey 3 age composition" << endl << srvp3 << endl;
  report << "Survey 3 expected age composition" << endl << Esrvp3 << endl;
  report << "Survey 3 observed and expected age comp" << endl << res_srv3 << endl;
  report << "Survey 3 Pearson residuals age comp" << endl << pearson_srv3 << endl; 
  report << "Survey 3 input N" << endl << multN_srv3 << endl;
  report << "Survey 3 effective N age comp" << endl << effN_srv3 << endl;   
  report << "Survey 3 length composition likelihoods" << endl << llsrvlenp3 << endl;
  report << "Survey 3 length composition" << endl << srvlenp3 << endl;
  report << "Survey 3 expected length composition" << endl << Esrvlenp3 << endl;
  report << "Survey 3 observed and expected length comp" << endl << res_srv3len << endl;
  report << "Survey 3 Pearson residuals length comp" << endl << pearson_srv3len << endl << endl;
  report << "Survey 4 q" << endl << q4 << endl;
  report << "Survey 4 power" << endl << q4_pow << endl; 
  report << "Survey 4 expected index" << endl << Eindxsurv4 << endl;
  report << "Survey 4 RMSE" << endl << RMSE_srv4 << endl << endl;
  report << "Survey 5 q" << endl << q5 << endl;
  report << "Survey 5 power" << endl << q5_pow << endl; 
  report << "Survey 5 expected index" << endl << Eindxsurv5 << endl;
  report << "Survey 5 RMSE" << endl << RMSE_srv5 << endl << endl;
  report << "Survey 6 q" << endl << q6 << endl;
  report << "Survey 6 selectivity parameters" << endl;  
  report << log_slp1_srv6 << " " << inf1_srv6 << " " << log_slp2_srv6 << " " << inf2_srv6 << endl; 
  report << "Survey 6 selectivity" << endl << slctsrv6 << endl;
  report << "Survey 6 expected index" << endl << Eindxsurv6 << endl;
  report << "Survey 6 RMSE" << endl << RMSE_srv6 << endl;
  report << "Survey 6 age composition likelihoods" << endl << llsrvp6 << endl;
  report << "Survey 6 age composition" << endl << srvp6 << endl;
  report << "Survey 6 expected age composition" << endl << Esrvp6 << endl;
  report << "Survey 6 observed and expected age comp" << endl << res_srv6 << endl;
  report << "Survey 6 Pearson residuals age comp" << endl << pearson_srv6 << endl;  
  report << "Survey 6 input N" << endl << multN_srv6 << endl;  
  report << "Survey 6 effective N age comp" << endl << effN_srv6 << endl;   
  report << "Survey 6 length composition likelihoods" << endl << llsrvlenp6 << endl;
  report << "Survey 6 length composition" << endl << srvlenp6 << endl;
  report << "Survey 6 expected length composition" << endl << Esrvlenp6 << endl << endl;
  report << "Projection selectivity" << endl << slctfsh_proj << endl;
  report << "Projection recruits" << endl << recruit_proj << endl;
  report << "Projection log mean recruitment" << endl << log_recr_proj << endl;
  report << "Projection recruitment variability" << endl <<  sigmasq_recr<< endl;
  report << "Projection numbers at age" << endl << N_proj << endl;
  report << "Projection catch at age" << endl << C_proj << endl; 
  report << "Projection survey numbers at age" << endl << Nsrv_proj << endl;
  report << "Projection population weight at age" << endl <<  wt_pop_proj << endl;  
  report << "Projection spawning weight at age" << endl << wt_spawn_proj << endl;  
  report << "Projection fishery weight at age" << endl << wt_fsh_proj << endl; 
  report << "Projection fishery selectivity" << endl << slctfsh_proj << endl; 
  report << "Projection total catches" << endl << Ecattot_proj << endl;  
  report << "Projection summary biomass" << endl << Esumbio_proj << endl;  
  report << "Projection spawning biomass" << endl << Espawnbio_proj << endl;  
  report << "Projection survey biomass" << endl << Esrv_proj << endl;  
  report << "Projection Ftarget" << endl << Ftarget << endl;
  report << "Projection B40" << endl <<  B40 << endl;  
  report << "Projection fishing mortality" << endl <<  F_proj << endl;  
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_report1;
  pad_report1 = NULL;
}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
 arrmblsize = 3000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
