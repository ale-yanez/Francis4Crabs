#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc2.csv");
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <LAM.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  ntime.allocate("ntime");
  nedades.allocate("nedades");
  ntallas.allocate("ntallas");
  data.allocate(1,ntime,1,13,"data");
  edades.allocate(1,nedades,"edades");
  Tallas.allocate(1,ntallas,"Tallas");
  Ctot.allocate(1,4,1,ntime,1,ntallas,"Ctot");
  msex.allocate(1,ntallas,"msex");
  Wmed.allocate(1,2,1,ntallas,"Wmed");
  cvar.allocate(1,3,"cvar");
  dt.allocate(1,3,"dt");
  Par_bio.allocate(1,2,1,5,"Par_bio");
  cv_priors.allocate(1,5,"cv_priors");
 log_Lopriorm = log(Par_bio(1,3));
 log_cva_priorm = log(Par_bio(1,4));
 log_M_priorm = log(Par_bio(1,5));
 log_Lopriorh = log(Par_bio(2,3));
 log_cva_priorh = log(Par_bio(2,4));
 log_M_priorh = log(Par_bio(2,5));
  h.allocate("h");
  qcru.allocate("qcru");
  cv_qcru.allocate("cv_qcru");
 log_qc_prior = log(qcru);
  parS.allocate(1,4,1,2,"parS");
 log_A50fpriorm = log(parS(1,1));
 log_s1priorm = log(parS(1,2));
 log_A50fpriorh = log(parS(2,1));
 log_s1priorh = log(parS(2,2));
 log_A50cpriorh   = log(parS(3,1));
 log_s1priorch = log(parS(3,2));
 log_A50cpriorm = log(parS(4,1));
 log_s1priorcm = log(parS(4,2));
  nbloq_selflo.allocate("nbloq_selflo");
  ybloq_selflo.allocate(1,nbloq_selflo,"ybloq_selflo");
  nbloq_selcru.allocate("nbloq_selcru");
  ybloq_selcru.allocate(1,nbloq_selcru,"ybloq_selcru");
  nbloq_qflo.allocate("nbloq_qflo");
  ybloq_qflo.allocate(1,nbloq_qflo,"ybloq_qflo");
  nbloq_qcru.allocate("nbloq_qcru");
  ybloq_qcru.allocate(1,nbloq_qcru,"ybloq_qcru");
  opt_qf.allocate("opt_qf");
  opt_qc.allocate("opt_qc");
  optSf_fase.allocate("optSf_fase");
  optSc_fase.allocate("optSc_fase");
  opt_Lo.allocate("opt_Lo");
  opt_cva.allocate("opt_cva");
  opt_M.allocate("opt_M");
  opt_F.allocate("opt_F");
  opt_devRt.allocate("opt_devRt");
  opt_devNo.allocate("opt_devNo");
  opt_Rm.allocate("opt_Rm");
  opt_Fpbr.allocate("opt_Fpbr");
  npbr.allocate("npbr");
  pbr.allocate(1,npbr,"pbr");
  ntime_sim.allocate("ntime_sim");
  pB.allocate("pB");
  opt_sim.allocate("opt_sim");
}

void model_parameters::initializationfunction(void)
{
  log_Lom.set_initial_value(log_Lopriorm);
  log_Loh.set_initial_value(log_Lopriorh);
  log_cv_edadm.set_initial_value(log_cva_priorm);
  log_cv_edadh.set_initial_value(log_cva_priorh);
  log_pRm.set_initial_value(-0.69314);
  log_A50flom.set_initial_value(log_A50fpriorm);
  log_A50floh.set_initial_value(log_A50fpriorh);
  log_sdA50flom.set_initial_value(log_s1priorm);
  log_sdA50floh.set_initial_value(log_s1priorh);
  log_A50crum.set_initial_value(log_A50cpriorm);
  log_A50cruh.set_initial_value(log_A50cpriorh);
  log_sdA50crum.set_initial_value(log_s1priorcm);
  log_sdA50cruh.set_initial_value(log_s1priorch);
  log_Mm.set_initial_value(log_M_priorm);
  log_Mh.set_initial_value(log_M_priorh);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_A50flom.allocate(1,nbloq_selflo,0.67,1.94,optSf_fase,"log_A50flom");
  log_sdA50flom.allocate(1,nbloq_selflo,optSf_fase,"log_sdA50flom");
  log_A50floh.allocate(1,nbloq_selflo,0.67,1.94,optSf_fase,"log_A50floh");
  log_sdA50floh.allocate(1,nbloq_selflo,optSf_fase,"log_sdA50floh");
  log_A50cruh.allocate(1,nbloq_selcru,0.67,1.94,optSc_fase,"log_A50cruh");
  log_sdA50cruh.allocate(1,nbloq_selcru,optSc_fase,"log_sdA50cruh");
  log_A50crum.allocate(1,nbloq_selcru,0.67,1.94,optSc_fase,"log_A50crum");
  log_sdA50crum.allocate(1,nbloq_selcru,optSc_fase,"log_sdA50crum");
  log_Ro.allocate(1,"log_Ro");
  log_pRm.allocate(-2.3,-0.1,opt_Rm,"log_pRm");
  dev_log_Ro.allocate(1,ntime,-10,10,opt_devRt,"dev_log_Ro");
  dev_log_Nom.allocate(1,nedades,-10,10,opt_devNo,"dev_log_Nom");
  dev_log_Noh.allocate(1,nedades,-10,10,opt_devNo,"dev_log_Noh");
  log_Fm.allocate(1,ntime,-20,-0.2,opt_F,"log_Fm");
  log_Fh.allocate(1,ntime,-20,-0.2,opt_F,"log_Fh");
  log_qflo.allocate(1,nbloq_qflo,opt_qf,"log_qflo");
  log_qcru.allocate(1,nbloq_qcru,opt_qc,"log_qcru");
  log_Lom.allocate(opt_Lo,"log_Lom");
  log_cv_edadm.allocate(opt_cva,"log_cv_edadm");
  log_Loh.allocate(opt_Lo,"log_Loh");
  log_cv_edadh.allocate(opt_cva,"log_cv_edadh");
  log_Mh.allocate(opt_M,"log_Mh");
  log_Mm.allocate(opt_M,"log_Mm");
  log_Fref.allocate(1,npbr,opt_Fpbr,"log_Fref");
  BMflo.allocate(1,ntime,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  BMcru.allocate(1,ntime,"BMcru");
  #ifndef NO_AD_INITIALIZE
    BMcru.initialize();
  #endif
  Brec.allocate(1,ntime,"Brec");
  #ifndef NO_AD_INITIALIZE
    Brec.initialize();
  #endif
  likeval.allocate(1,20,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  Neqm.allocate(1,nedades,"Neqm");
  #ifndef NO_AD_INITIALIZE
    Neqm.initialize();
  #endif
  Neqh.allocate(1,nedades,"Neqh");
  #ifndef NO_AD_INITIALIZE
    Neqh.initialize();
  #endif
  Rpred.allocate(1,ntime,"Rpred");
  #ifndef NO_AD_INITIALIZE
    Rpred.initialize();
  #endif
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_anos.allocate(1,ntime,"Unos_anos");
  #ifndef NO_AD_INITIALIZE
    Unos_anos.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  mu_edadm.allocate(1,nedades,"mu_edadm");
  #ifndef NO_AD_INITIALIZE
    mu_edadm.initialize();
  #endif
  mu_edadh.allocate(1,nedades,"mu_edadh");
  #ifndef NO_AD_INITIALIZE
    mu_edadh.initialize();
  #endif
  sigma_edadm.allocate(1,nedades,"sigma_edadm");
  #ifndef NO_AD_INITIALIZE
    sigma_edadm.initialize();
  #endif
  sigma_edadh.allocate(1,nedades,"sigma_edadh");
  #ifndef NO_AD_INITIALIZE
    sigma_edadh.initialize();
  #endif
  BDo.allocate(1,ntime,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  No.allocate(1,nedades,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  prior.allocate(1,7,"prior");
  #ifndef NO_AD_INITIALIZE
    prior.initialize();
  #endif
  prop_hpred.allocate(1,ntime,"prop_hpred");
  #ifndef NO_AD_INITIALIZE
    prop_hpred.initialize();
  #endif
  yrs.allocate(1,ntime,"yrs");
  #ifndef NO_AD_INITIALIZE
    yrs.initialize();
  #endif
  Desemb.allocate(1,ntime,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  CPUE.allocate(1,ntime,"CPUE");
  #ifndef NO_AD_INITIALIZE
    CPUE.initialize();
  #endif
  Bcru.allocate(1,ntime,"Bcru");
  #ifndef NO_AD_INITIALIZE
    Bcru.initialize();
  #endif
  prop_h.allocate(1,ntime,"prop_h");
  #ifndef NO_AD_INITIALIZE
    prop_h.initialize();
  #endif
  Lobs.allocate(1,ntime,"Lobs");
  #ifndef NO_AD_INITIALIZE
    Lobs.initialize();
  #endif
  Lpred.allocate(1,ntime,"Lpred");
  #ifndef NO_AD_INITIALIZE
    Lpred.initialize();
  #endif
  cv_index.allocate(1,4,1,ntime,"cv_index");
  #ifndef NO_AD_INITIALIZE
    cv_index.initialize();
  #endif
  nm_sex.allocate(1,4,1,ntime,"nm_sex");
  #ifndef NO_AD_INITIALIZE
    nm_sex.initialize();
  #endif
  Sel_crum.allocate(1,ntime,1,nedades,"Sel_crum");
  #ifndef NO_AD_INITIALIZE
    Sel_crum.initialize();
  #endif
  Sel_cruh.allocate(1,ntime,1,nedades,"Sel_cruh");
  #ifndef NO_AD_INITIALIZE
    Sel_cruh.initialize();
  #endif
  S1.allocate(1,nbloq_selflo,1,nedades,"S1");
  #ifndef NO_AD_INITIALIZE
    S1.initialize();
  #endif
  S2.allocate(1,nbloq_selflo,1,nedades,"S2");
  #ifndef NO_AD_INITIALIZE
    S2.initialize();
  #endif
  S3.allocate(1,nbloq_selcru,1,nedades,"S3");
  #ifndef NO_AD_INITIALIZE
    S3.initialize();
  #endif
  S4.allocate(1,nbloq_selcru,1,nedades,"S4");
  #ifndef NO_AD_INITIALIZE
    S4.initialize();
  #endif
  Sel_m.allocate(1,ntime,1,nedades,"Sel_m");
  #ifndef NO_AD_INITIALIZE
    Sel_m.initialize();
  #endif
  Sel_h.allocate(1,ntime,1,nedades,"Sel_h");
  #ifndef NO_AD_INITIALIZE
    Sel_h.initialize();
  #endif
  Fm.allocate(1,ntime,1,nedades,"Fm");
  #ifndef NO_AD_INITIALIZE
    Fm.initialize();
  #endif
  Fh.allocate(1,ntime,1,nedades,"Fh");
  #ifndef NO_AD_INITIALIZE
    Fh.initialize();
  #endif
  Zm.allocate(1,ntime,1,nedades,"Zm");
  #ifndef NO_AD_INITIALIZE
    Zm.initialize();
  #endif
  Zh.allocate(1,ntime,1,nedades,"Zh");
  #ifndef NO_AD_INITIALIZE
    Zh.initialize();
  #endif
  Sm.allocate(1,ntime,1,nedades,"Sm");
  #ifndef NO_AD_INITIALIZE
    Sm.initialize();
  #endif
  Sh.allocate(1,ntime,1,nedades,"Sh");
  #ifndef NO_AD_INITIALIZE
    Sh.initialize();
  #endif
  Nm.allocate(1,ntime,1,nedades,"Nm");
  #ifndef NO_AD_INITIALIZE
    Nm.initialize();
  #endif
  Nh.allocate(1,ntime,1,nedades,"Nh");
  #ifndef NO_AD_INITIALIZE
    Nh.initialize();
  #endif
  NM.allocate(1,ntime,1,nedades,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  NMD.allocate(1,ntime,1,ntallas,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  NDv.allocate(1,ntime,1,ntallas,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  Nrec.allocate(1,ntime,1,ntallas,"Nrec");
  #ifndef NO_AD_INITIALIZE
    Nrec.initialize();
  #endif
  NVflo_m.allocate(1,ntime,1,ntallas,"NVflo_m");
  #ifndef NO_AD_INITIALIZE
    NVflo_m.initialize();
  #endif
  NVflo_h.allocate(1,ntime,1,ntallas,"NVflo_h");
  #ifndef NO_AD_INITIALIZE
    NVflo_h.initialize();
  #endif
  NVcru_m.allocate(1,ntime,1,ntallas,"NVcru_m");
  #ifndef NO_AD_INITIALIZE
    NVcru_m.initialize();
  #endif
  NVcru_h.allocate(1,ntime,1,ntallas,"NVcru_h");
  #ifndef NO_AD_INITIALIZE
    NVcru_h.initialize();
  #endif
  pred_Ctotm.allocate(1,ntime,1,ntallas,"pred_Ctotm");
  #ifndef NO_AD_INITIALIZE
    pred_Ctotm.initialize();
  #endif
  pred_Ctot_am.allocate(1,ntime,1,nedades,"pred_Ctot_am");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot_am.initialize();
  #endif
  pred_Ctoth.allocate(1,ntime,1,ntallas,"pred_Ctoth");
  #ifndef NO_AD_INITIALIZE
    pred_Ctoth.initialize();
  #endif
  pred_Ctot_ah.allocate(1,ntime,1,nedades,"pred_Ctot_ah");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot_ah.initialize();
  #endif
  pobs_m.allocate(1,ntime,1,ntallas,"pobs_m");
  #ifndef NO_AD_INITIALIZE
    pobs_m.initialize();
  #endif
  ppred_m.allocate(1,ntime,1,ntallas,"ppred_m");
  #ifndef NO_AD_INITIALIZE
    ppred_m.initialize();
  #endif
  pobs_h.allocate(1,ntime,1,ntallas,"pobs_h");
  #ifndef NO_AD_INITIALIZE
    pobs_h.initialize();
  #endif
  ppred_h.allocate(1,ntime,1,ntallas,"ppred_h");
  #ifndef NO_AD_INITIALIZE
    ppred_h.initialize();
  #endif
  pobsc_m.allocate(1,ntime,1,ntallas,"pobsc_m");
  #ifndef NO_AD_INITIALIZE
    pobsc_m.initialize();
  #endif
  ppredc_m.allocate(1,ntime,1,ntallas,"ppredc_m");
  #ifndef NO_AD_INITIALIZE
    ppredc_m.initialize();
  #endif
  pobsc_h.allocate(1,ntime,1,ntallas,"pobsc_h");
  #ifndef NO_AD_INITIALIZE
    pobsc_h.initialize();
  #endif
  ppredc_h.allocate(1,ntime,1,ntallas,"ppredc_h");
  #ifndef NO_AD_INITIALIZE
    ppredc_h.initialize();
  #endif
  Prob_talla_m.allocate(1,nedades,1,ntallas,"Prob_talla_m");
  #ifndef NO_AD_INITIALIZE
    Prob_talla_m.initialize();
  #endif
  Prob_talla_h.allocate(1,nedades,1,ntallas,"Prob_talla_h");
  #ifndef NO_AD_INITIALIZE
    Prob_talla_h.initialize();
  #endif
  P1.allocate(1,nedades,1,ntallas,"P1");
  #ifndef NO_AD_INITIALIZE
    P1.initialize();
  #endif
  P2.allocate(1,nedades,1,ntallas,"P2");
  #ifndef NO_AD_INITIALIZE
    P2.initialize();
  #endif
  P3.allocate(1,nedades,1,ntallas,"P3");
  #ifndef NO_AD_INITIALIZE
    P3.initialize();
  #endif
  Nv.allocate(1,ntime,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  NMDv.allocate(1,ntime,1,nedades,"NMDv");
  #ifndef NO_AD_INITIALIZE
    NMDv.initialize();
  #endif
  suma1.allocate("suma1");
  #ifndef NO_AD_INITIALIZE
  suma1.initialize();
  #endif
  suma2.allocate("suma2");
  #ifndef NO_AD_INITIALIZE
  suma2.initialize();
  #endif
  suma3.allocate("suma3");
  #ifndef NO_AD_INITIALIZE
  suma3.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  suma5.allocate("suma5");
  #ifndef NO_AD_INITIALIZE
  suma5.initialize();
  #endif
  suma6.allocate("suma6");
  #ifndef NO_AD_INITIALIZE
  suma6.initialize();
  #endif
  suma7.allocate("suma7");
  #ifndef NO_AD_INITIALIZE
  suma7.initialize();
  #endif
  penalty.allocate("penalty");
  #ifndef NO_AD_INITIALIZE
  penalty.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Linfm.allocate("Linfm");
  #ifndef NO_AD_INITIALIZE
  Linfm.initialize();
  #endif
  km.allocate("km");
  #ifndef NO_AD_INITIALIZE
  km.initialize();
  #endif
  Linfh.allocate("Linfh");
  #ifndef NO_AD_INITIALIZE
  Linfh.initialize();
  #endif
  kh.allocate("kh");
  #ifndef NO_AD_INITIALIZE
  kh.initialize();
  #endif
  Mm.allocate("Mm");
  #ifndef NO_AD_INITIALIZE
  Mm.initialize();
  #endif
  Mh.allocate("Mh");
  #ifndef NO_AD_INITIALIZE
  Mh.initialize();
  #endif
  BDp.allocate("BDp");
  #ifndef NO_AD_INITIALIZE
  BDp.initialize();
  #endif
  Npplus.allocate("Npplus");
  #ifndef NO_AD_INITIALIZE
  Npplus.initialize();
  #endif
  Bph.allocate("Bph");
  #ifndef NO_AD_INITIALIZE
  Bph.initialize();
  #endif
  Bpm.allocate("Bpm");
  #ifndef NO_AD_INITIALIZE
  Bpm.initialize();
  #endif
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  nm2.allocate("nm2");
  #ifndef NO_AD_INITIALIZE
  nm2.initialize();
  #endif
  cuenta2.allocate("cuenta2");
  #ifndef NO_AD_INITIALIZE
  cuenta2.initialize();
  #endif
  nm3.allocate("nm3");
  #ifndef NO_AD_INITIALIZE
  nm3.initialize();
  #endif
  cuenta3.allocate("cuenta3");
  #ifndef NO_AD_INITIALIZE
  cuenta3.initialize();
  #endif
  nm4.allocate("nm4");
  #ifndef NO_AD_INITIALIZE
  nm4.initialize();
  #endif
  cuenta4.allocate("cuenta4");
  #ifndef NO_AD_INITIALIZE
  cuenta4.initialize();
  #endif
  Nph.allocate(1,nedades,"Nph");
  #ifndef NO_AD_INITIALIZE
    Nph.initialize();
  #endif
  Zpbrh.allocate(1,nedades,"Zpbrh");
  #ifndef NO_AD_INITIALIZE
    Zpbrh.initialize();
  #endif
  Fpbrh.allocate(1,nedades,"Fpbrh");
  #ifndef NO_AD_INITIALIZE
    Fpbrh.initialize();
  #endif
  Sph.allocate(1,nedades,"Sph");
  #ifndef NO_AD_INITIALIZE
    Sph.initialize();
  #endif
  Npm.allocate(1,nedades,"Npm");
  #ifndef NO_AD_INITIALIZE
    Npm.initialize();
  #endif
  Zpbrm.allocate(1,nedades,"Zpbrm");
  #ifndef NO_AD_INITIALIZE
    Zpbrm.initialize();
  #endif
  Fpbrm.allocate(1,nedades,"Fpbrm");
  #ifndef NO_AD_INITIALIZE
    Fpbrm.initialize();
  #endif
  Spm.allocate(1,nedades,"Spm");
  #ifndef NO_AD_INITIALIZE
    Spm.initialize();
  #endif
  CTP.allocate(1,ntallas,"CTP");
  #ifndef NO_AD_INITIALIZE
    CTP.initialize();
  #endif
  NMDp.allocate(1,ntallas,"NMDp");
  #ifndef NO_AD_INITIALIZE
    NMDp.initialize();
  #endif
  YTP.allocate(1,ntime_sim,1,npbr,"YTP");
  #ifndef NO_AD_INITIALIZE
    YTP.initialize();
  #endif
  BTp.allocate(1,ntime_sim,1,npbr,"BTp");
  #ifndef NO_AD_INITIALIZE
    BTp.initialize();
  #endif
  BD_lp.allocate("BD_lp");
  #ifndef NO_AD_INITIALIZE
  BD_lp.initialize();
  #endif
  ratio_pbr.allocate(1,npbr,"ratio_pbr");
  #ifndef NO_AD_INITIALIZE
    ratio_pbr.initialize();
  #endif
  Nvp.allocate(1,nedades,"Nvp");
  #ifndef NO_AD_INITIALIZE
    Nvp.initialize();
  #endif
  Nvplus.allocate("Nvplus");
  #ifndef NO_AD_INITIALIZE
  Nvplus.initialize();
  #endif
  SDvp.allocate(1,ntime_sim,"SDvp");
  #ifndef NO_AD_INITIALIZE
    SDvp.initialize();
  #endif
  CBA.allocate(1,npbr,"CBA");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  pred_CPUE.allocate(1,ntime,"pred_CPUE");
  pred_Bcru.allocate(1,ntime,"pred_Bcru");
  pred_Desemb.allocate(1,ntime,"pred_Desemb");
  BD.allocate(1,ntime,"BD");
  BT.allocate(1,ntime,"BT");
  BV.allocate(1,ntime,"BV");
  RPR.allocate(1,ntime,"RPR");
  SSBo.allocate("SSBo");
  RPRp.allocate(1,npbr,"RPRp");
  Restim.allocate(1,ntime,"Restim");
  RPRlp.allocate(1,ntime,"RPRlp");
  SSBp.allocate(1,ntime_sim,1,npbr,"SSBp");
  Lmf_obs.allocate(1,ntime,"Lmf_obs");
  #ifndef NO_AD_INITIALIZE
    Lmf_obs.initialize();
  #endif
  Lmf_pred.allocate(1,ntime,"Lmf_pred");
  Lhf_obs.allocate(1,ntime,"Lhf_obs");
  #ifndef NO_AD_INITIALIZE
    Lhf_obs.initialize();
  #endif
  Lhf_pred.allocate(1,ntime,"Lhf_pred");
  Lmc_obs.allocate(1,ntime,"Lmc_obs");
  #ifndef NO_AD_INITIALIZE
    Lmc_obs.initialize();
  #endif
  Lmc_pred.allocate(1,ntime,"Lmc_pred");
  Lhc_obs.allocate(1,ntime,"Lhc_obs");
  #ifndef NO_AD_INITIALIZE
    Lhc_obs.initialize();
  #endif
  Lhc_pred.allocate(1,ntime,"Lhc_pred");
  Frpr.allocate(1,ntime,"Frpr");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 yrs=column(data,1);
 Desemb=column(data,2);
 CPUE=column(data,3);
 Bcru=column(data,4);
 prop_h=column(data,5);
 for (int i=1;i<=4;i++){
 cv_index(i)=column(data,i+5);
 nm_sex(i)=column(data,i+9);
 }
 Linfm=Par_bio(1,1);
 km=Par_bio(1,2);
 Linfh=Par_bio(2,1);
 kh=Par_bio(2,2);
 // Mm=Par_bio(1,5);
 // Mh=Par_bio(2,5);
 
 Unos_edad=1;// lo uso en  operaciones matriciales con la edad
 Unos_anos=1;// lo uso en operaciones matriciales con el año
 Unos_tallas=1;// lo uso en operaciones matriciales con el año
  reporte_mcmc=0;
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-1,1.e-01,1.e-03,1e-3,1e-4}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100,100,200,300,3500}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::userfunction(void)
{
  f =0.0;
 Eval_prob_talla_edad();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_deinteres();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_PBR();
 Eval_logverosim();
 Eval_funcion_objetivo();
 if(last_phase){
 Eval_CTP();
 Eval_mcmc();}
}

void model_parameters::Eval_prob_talla_edad(void)
{
 int i, j;
 mu_edadm(1)=exp(log_Lom);
 for (i=2;i<=nedades;i++)
  {
  mu_edadm(i)=Linfm*(1-exp(-km))+exp(-km)*mu_edadm(i-1);
  }
 sigma_edadm=exp(log_cv_edadm)*mu_edadm;
 mu_edadh(1)=exp(log_Loh);
 for (i=2;i<=nedades;i++)
  {
  mu_edadh(i)=Linfh*(1-exp(-kh))+exp(-kh)*mu_edadh(i-1);
  }
 sigma_edadh=exp(log_cv_edadh)*mu_edadh;
  Prob_talla_m = ALK( mu_edadm, sigma_edadm, Tallas);
  Prob_talla_h = ALK( mu_edadh, sigma_edadh, Tallas);
}

dvar_matrix model_parameters::ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
{
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
}

void model_parameters::Eval_selectividad(void)
{
 int i,j;
 for (j=1;j<=nbloq_selflo;j++){
 S1(j)=1/(1+exp(-log(19)*(edades-exp(log_A50flom(j)))/exp(log_sdA50flom(j))));//machos
 S2(j)=1/(1+exp(-log(19)*(edades-exp(log_A50floh(j)))/exp(log_sdA50floh(j))));//hembras
 }
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloq_selflo;j++){
              if (yrs(i)>=ybloq_selflo(j)){
                Sel_m(i)=S1(j);//machos
                Sel_h(i)=S2(j);} //hembras
       }
   }
 // CRUCEROS
 // por defecto los mismos que la flota
 //    Sel_crum=Sel_m;
 //   Sel_cruh=Sel_h;
    Sel_crum=1.0;
    Sel_cruh=1.0;
 if(active(log_A50crum)){
 for (j=1;j<=nbloq_selcru;j++){
 S3(j)=1/(1+exp(-log(19)*(edades-exp(log_A50crum(j)))/exp(log_sdA50crum(j))));
 S4(j)=1/(1+exp(-log(19)*(edades-exp(log_A50cruh(j)))/exp(log_sdA50cruh(j))));
 }
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloq_selcru;j++){
              if (yrs(i)>=ybloq_selcru(j)){
                Sel_crum(i)=S3(j);
                Sel_cruh(i)=S4(j);}
       }
   }
 }
  /*
 for (j=1;j<=nbloq_selflo;j++){
 S1(j)=exp(-0.5*square(Tallas-exp(log_A50flom(j)))/square(exp(log_sigma1m(j))));//machos
 S2(j)=exp(-0.5*square(Tallas-exp(log_A50floh(j)))/square(exp(log_sigma1h(j))));//hembras
    for (i=1;i<=ntallas;i++){
      if(Tallas(i)>=exp(log_A50flom(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_A50flom(j)))/square(exp(log_sigma2m(j))));
      }
      if(Tallas(i)>=exp(log_A50floh(j))){
      S2(j,i)= exp(-0.5*square(Tallas(i)-exp(log_A50floh(j)))/square(exp(log_sigma2h(j))));
      }
 }}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloq_selflo;j++){
              if (yrs(i)>=ybloq_selflo(j)){
                Sel_m(i)=Prob_talla_m*S1(j);//machos
                Sel_h(i)=Prob_talla_h*S2(j);} //hembras
       }
   }
 // CRUCEROS
 // por defecto los mismos que la flota
 //    Sel_crum=Sel_m;
 //   Sel_cruh=Sel_h;
    Sel_crum=1.0;
    Sel_cruh=1.0;
 if(active(log_A50crum)){
 for (j=1;j<=nbloq_selcru;j++){
 S3(j)=exp(-0.5*square(Tallas-exp(log_A50crum(j)))/square(exp(log_sigma1cm(j))));
 S4(j)=exp(-0.5*square(Tallas-exp(log_A50cruh(j)))/square(exp(log_sigma1ch(j))));
    for (i=1;i<=ntallas;i++){
      if(Tallas(i)>=exp(log_A50crum(j))){
      S3(j,i)= exp(-0.5*square(Tallas(i)-exp(log_A50crum(j)))/square(exp(log_sigma2cm(j))));
      }
      if(Tallas(i)>=exp(log_A50cruh(j))){
      S4(j,i)= exp(-0.5*square(Tallas(i)-exp(log_A50cruh(j)))/square(exp(log_sigma2ch(j))));
      }
 }}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloq_selcru;j++){
              if (yrs(i)>=ybloq_selcru(j)){
                Sel_crum(i)=Prob_talla_m*S3(j);
                Sel_cruh(i)=Prob_talla_h*S4(j);}
       }
   }
 }
 */
}

void model_parameters::Eval_mortalidades(void)
{
 Mm=exp(log_Mm);
 Mh=exp(log_Mh);
 Fm=elem_prod(Sel_m,outer_prod(mfexp(log_Fm),Unos_edad));
 Fh=elem_prod(Sel_h,outer_prod(mfexp(log_Fh),Unos_edad));
 Zm=Fm+Mm;
 Zh=Fh+Mh;
 Sm=mfexp(-1.0*Zm);
 Sh=mfexp(-1.0*Zh);
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
 // Biomasa desovante virgen de largo plazo
 No(1)=exp(log_Ro); // hembras
 for (int j=2;j<=nedades;j++)
 {
   No(j)=No(j-1)*exp(-1.*Mh);
 }
   No(nedades)=No(nedades)*exp(-1.*Mh)/(1-exp(-1.*Mh));
  SSBo=sum(elem_prod(No*exp(-dt(1)*Mh)*Prob_talla_h,elem_prod(msex,Wmed(2))));
 alfa=4*h*exp(log_Ro)/(5*h-1);//
 beta=(1-h)*SSBo/(5*h-1);//  
 Neqh(1)=mfexp(log_Ro);//hembras
 for (j=2;j<=nedades;j++)
 {
   Neqh(j)=Neqh(j-1)*exp(-Zh(1,j-1));
 }
   Neqh(nedades)=Neqh(nedades)*exp(-1.*Zh(1,nedades))/(1-exp(-1.*Zh(1,nedades)));
 Neqm(1)=mfexp(log_Ro) * (exp(log_pRm)/(1-exp(log_pRm)));//machos
 for (j=2;j<=nedades;j++)
 {
   Neqm(j)=Neqm(j-1)*exp(-Zm(1,j-1));
 }
   //Neqm(nedades)+=Neqm(nedades)*exp(-1.*Zm(1,nedades));//
   Neqm(nedades)=Neqm(nedades)*exp(-1.*Zm(1,nedades))/(1-exp(-1.*Zm(1,nedades)));
 Nm(1)=elem_prod(Neqm,exp(dev_log_Nom));
 Nh(1)=elem_prod(Neqh,exp(dev_log_Noh));
 BD(1)=sum(elem_prod(elem_prod(Nh(1),exp(-dt(1)*Zh(1)))*Prob_talla_h,elem_prod(msex,Wmed(2))));
 Rpred(1)=mfexp(log_Ro);//
 for (i=1;i<ntime;i++)
 {
     Rpred(i+1)=mfexp(log_Ro);// 
     if(i>=edades(1)){
     Rpred(i+1)=alfa*BD(i-edades(1)+1)/(beta+BD(i-edades(1)+1));}// Reclutamiento estimado por un modelo B&H hembras
     Nm(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i))*exp(log_pRm)/(1-exp(log_pRm));  // Reclutas machos   
     Nh(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i));// Reclutas hembras
     Restim=column(Nh,1);
     Nm(i+1)(2,nedades)=++elem_prod(Nm(i)(1,nedades-1),Sm(i)(1,nedades-1));
     Nm(i+1,nedades)=Nm(i+1,nedades)+Nm(i,nedades)*Sm(i,nedades);// grupo plus
     Nh(i+1)(2,nedades)=++elem_prod(Nh(i)(1,nedades-1),Sh(i)(1,nedades-1));
     Nh(i+1,nedades)=Nh(i+1,nedades)+Nh(i,nedades)*Sh(i,nedades);// grupo plus
     BD(i+1)=sum(elem_prod(elem_prod(Nh(i+1),exp(-dt(1)*Zh(i+1)))*Prob_talla_h,elem_prod(msex,Wmed(2))));
 }
}

void model_parameters::Eval_deinteres(void)
{
 Nv=Nh;// solo para empezar los calculos
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*Mh);
     Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*Mh);// grupo plus
 }
 NDv=elem_prod((Nv*exp(-dt(1)*Mh))*Prob_talla_h,outer_prod(Unos_anos,msex));
 BDo=NDv*Wmed(2);
 RPR=elem_div(BD,BDo);
 RPRlp=BD/SSBo;
}

void model_parameters::Eval_biomasas(void)
{
 NMD=elem_prod(Nh,mfexp(-dt(1)*Zh))*Prob_talla_h;
 NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));
 NVflo_m=elem_prod(elem_prod(Nm,mfexp(-dt(2)*(Zm))),Sel_m)*Prob_talla_m;
 NVflo_h=elem_prod(elem_prod(Nh,mfexp(-dt(2)*(Zh))),Sel_h)*Prob_talla_h;
 NVcru_m=elem_prod(elem_prod(Nm,mfexp(-dt(3)*(Zm))),Sel_crum)*Prob_talla_m;
 NVcru_h=elem_prod(elem_prod(Nh,mfexp(-dt(3)*(Zh))),Sel_cruh)*Prob_talla_h;
 BD=NMD*Wmed(2);
 BMflo=NVflo_m*Wmed(1)+NVflo_h*Wmed(2);
 BMcru=NVcru_m*Wmed(1)+NVcru_h*Wmed(2);
 BV=BMflo;
 BT=(Nm*Prob_talla_m)*Wmed(1)+(Nh*Prob_talla_h)*Wmed(2);  
}

void model_parameters::Eval_capturas_predichas(void)
{
 pred_Ctot_am=elem_prod(elem_div(Fm,Zm),elem_prod(1.-Sm,Nm));
 pred_Ctotm=pred_Ctot_am*Prob_talla_m;
 pred_Ctot_ah=elem_prod(elem_div(Fh,Zh),elem_prod(1.-Sh,Nh));
 pred_Ctoth=pred_Ctot_ah*Prob_talla_h;
 prop_hpred = elem_div(rowsum(pred_Ctoth),rowsum(pred_Ctoth+pred_Ctotm+1e-10));
 pred_Desemb=pred_Ctotm*Wmed(1)+pred_Ctoth*Wmed(2);
 pobs_m=elem_div(Ctot(1),outer_prod(rowsum(Ctot(1)+1e-10),Unos_tallas));
 ppred_m=elem_div(pred_Ctotm,outer_prod(rowsum(pred_Ctotm+1e-10),Unos_tallas));
 pobs_h=elem_div(Ctot(2),outer_prod(rowsum(Ctot(2)+1e-10),Unos_tallas));
 ppred_h=elem_div(pred_Ctoth,outer_prod(rowsum(pred_Ctoth+1e-10),Unos_tallas));
 pobsc_m=elem_div(Ctot(3),outer_prod(rowsum(Ctot(3)+1e-10),Unos_tallas));
 ppredc_m=elem_div(NVcru_m,outer_prod(rowsum(NVcru_m+1e-10),Unos_tallas));
 pobsc_h=elem_div(Ctot(4),outer_prod(rowsum(Ctot(4)+1e-10),Unos_tallas));
 ppredc_h=elem_div(NVcru_h,outer_prod(rowsum(NVcru_h+1e-10),Unos_tallas));
 Lmf_obs  = Tallas*trans(pobs_m);
 Lmf_pred = Tallas*trans(ppred_m);
 Lhf_obs  = Tallas*trans(pobs_h);
 Lhf_pred = Tallas*trans(ppred_h);
 Lmc_obs  = Tallas*trans(pobsc_m);
 Lmc_pred = Tallas*trans(ppredc_m);
 Lhc_obs  = Tallas*trans(pobsc_h);
 Lhc_pred = Tallas*trans(ppredc_h);
}

void model_parameters::Eval_indices(void)
{
   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nbloq_qflo;j++){
              if (yrs(i)>=ybloq_qflo(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*BMflo(i);}
       }
   }
   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nbloq_qcru;j++){
              if (yrs(i)>=ybloq_qcru(j)){
                 pred_Bcru(i)=exp(log_qcru(j))*BMcru(i);}
       }
   }
}

void model_parameters::Eval_PBR(void)
{
 for (int i=1;i<=npbr;i++){
 Fpbrh=Sel_h(ntime)*exp(log_Fref(i));
 Zpbrh=Fpbrh+Mh;
 Neqh(1)=mfexp(log_Ro);//hembras
 for (int j=2;j<=nedades;j++)
 {
   Neqh(j)=Neqh(j-1)*exp(-Zpbrh(j-1));
 }
   Neqh(nedades)=Neqh(nedades)*exp(-1.*Zpbrh(nedades))/(1-exp(-1.*Zpbrh(nedades))); // MODIFICAR POR LA OTRA FORMA
 BD_lp=sum(elem_prod(elem_prod(Neqh,exp(-dt(1)*Zpbrh))*Prob_talla_h,elem_prod(msex,Wmed(2))));
 ratio_pbr(i)=BD_lp/SSBo;
 }
 Frpr=mfexp(log_Fh)/mfexp(log_Fref(3));
}

void model_parameters::Eval_logverosim(void)
{
 int i;
 suma1=0; suma2=0; suma3=0; penalty=0;
 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}
  if (Bcru(i)>0){
    suma2+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(3,i));}
  if (prop_h(i)>0){
    suma3+=square(log(prop_h(i)/prop_hpred(i))*1/cv_index(4,i));}
 }
}

void model_parameters::Eval_funcion_objetivo(void)
{
 suma4=0; suma5=0; suma6=0; suma7=0; penalty=0;
 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//Crucero
 likeval(3)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb
 likeval(4)=0.5*suma3;// prop p_hembras
 likeval(5)=-1.*sum(nm_sex(1)*elem_prod(pobs_m,log(ppred_m)));
 likeval(6)=-1.*sum(nm_sex(2)*elem_prod(pobs_h,log(ppred_h)));
 likeval(7)=-1.*sum(nm_sex(3)*elem_prod(pobsc_m,log(ppredc_m)));
 likeval(8)=-1.*sum(nm_sex(4)*elem_prod(pobsc_h,log(ppredc_h)));
 if(active(dev_log_Ro)){
 likeval(9)=1./(2*square(cvar(1)))*norm2(dev_log_Ro);}
 if(active(dev_log_Nom)){
 likeval(10)=1./(2*square(cvar(2)))*norm2(dev_log_Nom);
 likeval(11)=1./(2*square(cvar(2)))*norm2(dev_log_Noh);}
 if (active(log_pRm)){
 likeval(16)=0.5/square(cvar(3))*square(log_pRm+0.69315);}
 if (active(log_Lom)){
 likeval(17)=0.5*square((log_Lopriorm-log_Lom)/cv_priors(3));
 likeval(18)=0.5*square((log_Lopriorh-log_Loh)/cv_priors(3));}
 if (active(log_cv_edadm)){
 likeval(19)=0.5*square((log_cv_edadm-log_cva_priorm)/cv_priors(4));
 likeval(20)=0.5*square((log_cv_edadh-log_cva_priorh)/cv_priors(4));}
 if(active(log_Mh)){
 penalty+=100*(square(log_M_priorh-log_Mh)+square(log_M_priorm-log_Mm));}
 if(active(log_qcru)){
 penalty+=0.5*norm2((log_qcru-log_qc_prior)/cv_qcru);}
 if(active(log_Fref)){
 penalty+=1000*norm2(ratio_pbr-pbr);}
 f=opt_sim*(sum(likeval)+penalty);
}

void model_parameters::Eval_CTP(void)
{
 // se considera el Fpbr de hembras como el representativo factor limitante
 for (int j=1;j<=npbr;j++){ // son # PBR only!
 Nph=Nh(ntime);
 Npm=Nm(ntime);
 Sph=Sh(ntime);
 Spm=Sm(ntime);
 BDp=BD(ntime);
 Fpbrh=Fh(ntime);//
 Fpbrm=Fm(ntime);//
 Zpbrh=Zh(ntime);
 Zpbrm=Zm(ntime);
 for (int i=1;i<=ntime_sim;i++)
 {
 Bph=sum(elem_prod(Nph*Prob_talla_h,Wmed(2)));
 Bpm=sum(elem_prod(Npm*Prob_talla_m,Wmed(1)));
 NMDp=elem_prod(Nph,mfexp(-dt(1)*(Zpbrh)))*Prob_talla_h;
 BDp=sum(elem_prod(elem_prod(NMDp,msex),Wmed(2)));
 CTP=elem_prod(elem_prod(elem_div(Fpbrh,Zpbrh),elem_prod(Nph,(1-Sph)))*Prob_talla_h,Wmed(2));
 CTP+=elem_prod(elem_prod(elem_div(Fpbrm,Zpbrm),elem_prod(Npm,(1-Spm)))*Prob_talla_m,Wmed(1));
 YTP(i,j)=sum(CTP);
 SSBp(i,j)=BDp;
 BTp(i,j)=Bph+Bpm;
 // año siguiente
 Npplus=Nph(nedades)*Sph(nedades);
 Nph(2,nedades)=++elem_prod(Nph(1,nedades-1),Sph(1,nedades-1));
 Nph(nedades)+=Npplus;
 Nph(1)=pB*exp(log_Ro);
 Npplus=Npm(nedades)*Spm(nedades);
 Npm(2,nedades)=++elem_prod(Npm(1,nedades-1),Spm(1,nedades-1));
 Npm(nedades)+=Npplus;
 Npm(1)=exp(log_Ro)*exp(log_pRm)/(1-exp(log_pRm));//;
 // Se considera el mismo F de hembras en los machos
 Fpbrh=Sel_h(ntime)*exp(log_Fref(j));//
 Fpbrm=Sel_m(ntime)*exp(log_Fref(j));//
 Zpbrh=Fpbrh+Mh;
 Zpbrm=Fpbrm+Mm;
 Sph=exp(-1.*Zpbrh);
 Spm=exp(-1.*Zpbrm);
 }}
 CBA=YTP(2);// es para el año proyectado
 // Rutina para la estimación de RPR
 Nvp=Nv(ntime);// toma la ultima estimación
 for (int i=1;i<=ntime_sim;i++)
  {
     Nvplus=Nvp(nedades)*exp(-1.0*Mh);
     Nvp(2,nedades)=++Nvp(1,nedades-1)*exp(-1.0*Mh);
     Nvp(nedades)+=Nvplus;
     Nvp(1)= exp(log_Ro);
     SDvp(i)=sum(elem_prod(Nvp*Prob_talla_h,elem_prod(Wmed(2),msex)));
  }
 for (int i=1;i<=npbr;i++)
 {
  RPRp(i)=SSBp(ntime_sim,i)/SDvp(ntime_sim);//
 }
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
 report << "YRS" << endl;
 report << yrs << endl;
 report << "CPUE" << endl;
 report << CPUE << endl;
 report << pred_CPUE << endl;
 report << "BCRU" << endl;
 report << Bcru << endl;
 report << pred_Bcru << endl;
 report << "Desemb" << endl;
 report << Desemb << endl;
 report << pred_Desemb << endl;
 report << "BD" << endl;
 report << BD << endl;
 report << "BT" << endl;
 report << BT << endl;
 report << "BV" << endl;
 report << BMflo << endl;
 report << "Rech_pre_est" << endl;
 report << Rpred<< endl;
 report << column(Nh,1)<< endl;
 report << "Fm_Fh" << endl;
 report << rowsum(Fm)/nedades<< endl;
 report << rowsum(Fh)/nedades<< endl;
 report << "PP_h_obs_pre" << endl;
 report << prop_h << endl;
 report << prop_hpred << endl;
 report << "Lm_obs_pred" << endl;
 report << Tallas*trans(pobs_m)<< endl;
 report << Tallas*trans(ppred_m)<< endl;
 report << "Lh_obs_pred" << endl;
 report << Tallas*trans(pobs_h)<< endl;
 report << Tallas*trans(ppred_h)<< endl;
 report << "Lmc_obs_est" << endl;
 report << Tallas*trans(pobsc_m)<< endl;
 report << Tallas*trans(ppredc_m)<< endl;
 report << "Lhc_obs_est" << endl;
 report << Tallas*trans(pobsc_h)<< endl;
 report << Tallas*trans(ppredc_h)<< endl;
 report << "Sflom_age" << endl;
 report << Sel_m << endl;
 report << "Sfloh_age" <<endl;
 report << Sel_h << endl;
 report << "Scrum_age" << endl;
 report << Sel_crum << endl;
 report << "Scruh_age" << endl;
 report << Sel_cruh << endl;
 report << "CAPTURAS" << endl;
 report << "pobs_mflo" << endl;
 report << (pobs_m)<< endl;
 report << "ppred_mflo" << endl;
 report << (ppred_m)<< endl;
 report << "pobs_hflo" << endl;
 report << (pobs_h)<< endl;
 report << "Ppred_hflo" << endl;
 report << (ppred_h)<< endl;
 report << "CRUCERO" << endl;
 report << "pobs_mcru" << endl;
 report << (pobsc_m)<< endl;
 report << "ppred_mcru" << endl;
 report << (ppredc_m)<< endl;
 report << "pobs_hcru" << endl;
 report << (pobsc_h)<< endl;
 report << "ppred_hcru" << endl;
 report << (ppredc_h)<< endl;
 report << "Abundancia a la edad"<< endl;
 report << "Nm"<< endl;
 report << Nm << endl;
 report << "Nh" << endl;
 report << Nh << endl;
 report << "Captura a la edad" <<endl;
 report << "Capt_mage" <<endl;
 report << pred_Ctot_am <<endl;
 report << "Capt_hage" <<endl;
 report << pred_Ctot_ah <<endl;
 report << "F_age" <<endl;
 report << "Fm_age" <<endl;
 report << Fm << endl;
 report << "Fh_age" << endl;
 report << Fh << endl;
 report << "BDo" << endl;
 report << BDo << endl;
 report << "BDoLP" << endl;
 report << SSBo << endl;
 report << "RPR" << endl;
 report << RPR << endl;
 report << "RPRlp" << endl;
 report << RPRlp << endl;
 report << "Lo_h"  << endl;
 report << exp(log_Loh) << endl;
 report << "cv_h" << endl;
 report << exp(log_cv_edadh) << endl;
 report << "Lo_m" << endl;
 report << exp(log_Lom) << endl;
 report << "cv_m" << endl;
 report << exp(log_cv_edadm) << endl;
 report << "Mm" << endl;
 report << Mm << endl;
 report << "Mh" << endl;
 report << Mh <<endl;
 report << "mu_edadm" << endl;
 report << mu_edadm << endl;
 report << "mu_edadh" << endl;
 report << mu_edadh << endl;
 report << "Prob_talla_m" << endl;
 report << Prob_talla_m << endl;
 report << "Prob_talla_h" << endl;
 report << Prob_talla_h << endl;
 report << "LIKE" <<endl; //CPUE, Crucero,Desemb, prop,prop_mflo, prop_hflo,pobsc_m, pobsc_h, Ro,No_m, No_h,Lo_m,Lo_h,cvage_m, cvage_h
 report << likeval << endl;
 report << "q_cru" <<endl;
 report << exp(log_qcru) << endl;
 report << "q_flo" <<endl;
 report << exp(log_qflo) << endl;
 report << "alfa" << endl;
 report << alfa << endl;
 report << "beta" <<endl;
 report << beta << endl;
 report << "ratio_obj" <<endl;
 report << ratio_pbr << endl;
 report << " Fpbr" << endl;
 report << exp(log_Fref) << endl;
 report << " RPR_lp" << endl;
 report << RPRp << endl;
 report << "BT_proy" << endl;
 report << BTp << endl;
 report << "BD_proy" << endl;
 report << SSBp << endl;
 report << "C_proy" << endl;
 report << YTP << endl;
 suma1=0; suma2=0;nm1=1;cuenta1=0;
  for (int i=1;i<=ntime;i++){ //
   if (sum(pobs_m(i))>0){
      suma1=sum(elem_prod(ppred_m(i),1-ppred_m(i)));
      suma2=norm2(pobs_m(i)-ppred_m(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}
 suma1=0;suma2=0;nm2=1;cuenta2=0;
  for (int i=1;i<=ntime;i++){ //
   if (sum(pobs_h(i))>0){
      suma1=sum(elem_prod(ppred_h(i),1-ppred_h(i)));
      suma2=norm2(pobs_h(i)-ppred_h(i));
      nm2=nm2*suma1/suma2;
      cuenta2+=1;
   }}
 suma1=0;suma2=0;nm3=1;cuenta3=0;
  for (int i=1;i<=ntime;i++){ //
   if (sum(pobsc_m(i))>0){
      suma1=sum(elem_prod(ppredc_m(i),1-ppredc_m(i)));
      suma2=norm2(pobsc_m(i)-ppredc_m(i));
      nm3=nm3*suma1/suma2;
      cuenta3+=1;
   }}
 suma1=0;suma2=0;nm4=1;cuenta4=0;
  for (int i=1;i<=ntime;i++){ //
   if (sum(pobs_m(i))>0){
      suma1=sum(elem_prod(ppredc_h(i),1-ppredc_h(i)));
      suma2=norm2(pobsc_h(i)-ppredc_h(i));
      nm4=nm4*suma1/suma2;
      cuenta4+=1;
   }}
 report << "Tamanho muestra ideal" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;
 report <<pow(nm2,1/cuenta2)<< endl;
 report <<pow(nm3,1/cuenta3)<< endl;
 report <<pow(nm4,1/cuenta4)<< endl;
}

void model_parameters::Eval_mcmc(void)
{
  if(reporte_mcmc == 0)
  mcmc_report<<"Bcru_last CTP1 CTP2 CTP3 CTP4 BDp1/BDlast BDp2/BDlast BDp3/BDlast BDp4/BDlast "<<endl;
  mcmc_report<<pred_Bcru(ntime)<<" "<<YTP(2,1)<<" "<<YTP(2,2)<<" "<<YTP(2,3)<<" "<<YTP(2,4)<<
     " "<<SSBp(ntime_sim,1)/BD(ntime)<<" "<<SSBp(ntime_sim,2)/BD(ntime)<<" "<<SSBp(ntime_sim,3)/BD(ntime)<<
     " "<<SSBp(ntime_sim,4)/BD(ntime)<<endl;
  reporte_mcmc++;
}

void model_parameters::final_calcs()
{
 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

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
 time(&start);
 arrmblsize = 50000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
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
