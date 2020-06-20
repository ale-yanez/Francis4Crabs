GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc2.csv");

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 50000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 

DATA_SECTION
 init_int ntime
 init_int nedades
 init_int ntallas

 init_matrix data(1,ntime,1,13);
  //Yrs/desem/cpue/Bcru/ph/cv_desem/cv_CPUE/cv_Cru/cv_ph/n_mf/n_hf/n_mc/n_hc
 init_vector edades(1,nedades);
 init_vector Tallas(1,ntallas);

 init_3darray Ctot(1,4,1,ntime,1,ntallas);

 init_vector msex(1,ntallas);
 init_matrix Wmed(1,2,1,ntallas);

//!! ad_comm::change_datafile_name("17.ctl");
 init_vector cvar(1,3);//# Coeficiente de variación de los desvios Rt, No y prop_machos en el Reclutamiento
 init_vector dt(1,3);
 init_matrix Par_bio(1,2,1,5); //#Loo, k, Lt(a=1), cv(edad), M; filas machos,hembras y prior si se estiman
 init_vector cv_priors(1,5);

  number log_Lopriorm
  number log_Lopriorh

  number log_cva_priorm
  number log_cva_priorh

  number log_M_priorm
  number log_M_priorh


  !! log_Lopriorm = log(Par_bio(1,3));
  !! log_cva_priorm = log(Par_bio(1,4));
  !! log_M_priorm = log(Par_bio(1,5));
	
  !! log_Lopriorh = log(Par_bio(2,3));
  !! log_cva_priorh = log(Par_bio(2,4));
  !! log_M_priorh = log(Par_bio(2,5));


 init_number h
 
 init_number qcru
 init_number cv_qcru

  number log_qc_prior
  !! log_qc_prior = log(qcru);

 init_matrix parS(1,4,1,2);//Selectividad flota y crucero (valores de partida y rango)
 
 number log_A50fpriorm //Selectividad flota
 number log_s1priorm
 number log_A50fpriorh
 number log_s1priorh

 !! log_A50fpriorm = log(parS(1,1));
 !! log_s1priorm = log(parS(1,2));
 !! log_A50fpriorh = log(parS(2,1));
 !! log_s1priorh = log(parS(2,2));

 number log_A50cpriorm
 number log_s1priorcm
 number log_A50cpriorh  //Selectividad cuceros
 number log_s1priorch

 !! log_A50cpriorh   = log(parS(3,1));
 !! log_s1priorch = log(parS(3,2));
 !! log_A50cpriorm = log(parS(4,1));
 !! log_s1priorcm = log(parS(4,2));

 init_int    nbloq_selflo // bloques selectividad flota
 init_vector ybloq_selflo(1,nbloq_selflo);

 init_int    nbloq_selcru // bloques selectividad crucero
 init_vector ybloq_selcru(1,nbloq_selcru);

 init_int    nbloq_qflo  //// bloques capturabilidad cpue
 init_vector ybloq_qflo(1,nbloq_qflo);

 init_int    nbloq_qcru   // bloques capturabilidad crucero
 init_vector ybloq_qcru(1,nbloq_qcru);

// Fases de estimacion
//_________________________________________________________
 init_int    opt_qf  // q flota
 init_int    opt_qc  // q crucero

 init_int    optSf_fase //fase sel flo
 init_int    optSc_fase //fase sel cru

 init_int    opt_Lo
 init_int    opt_cva
 init_int    opt_M

 init_int    opt_F
 init_int    opt_devRt
 init_int    opt_devNo
 init_int    opt_Rm
 init_int    opt_Fpbr


 init_int    npbr
 init_vector pbr(1,npbr);
 
 init_int ntime_sim

 init_number pB // Proyección de capturas ante distintos niveles de Reclutamiento (1.0 proporcional al reclutamiento medio)
 init_number opt_sim // Opción para simular o estimar(0=simula, 1=estima)
 int reporte_mcmc


INITIALIZATION_SECTION

  log_Lom        log_Lopriorm
  log_Loh        log_Lopriorh
	
  log_cv_edadm   log_cva_priorm // Parbio = 0.1
  log_cv_edadh   log_cva_priorh // Parbio = 0.1
  
	log_pRm        -0.69314

  log_A50flom        log_A50fpriorm 
  log_A50floh        log_A50fpriorh 
	  
  log_sdA50flom     log_s1priorm 
	log_sdA50floh     log_s1priorh 
	
  log_A50crum        log_A50cpriorm 
  log_A50cruh        log_A50cpriorh  
	
  log_sdA50crum     log_s1priorcm 
  log_sdA50cruh     log_s1priorch 

  log_Mm           log_M_priorm
  log_Mh           log_M_priorh


PARAMETER_SECTION

// init_bounded_vector log_L50f(1,nbloq_selflo,-5,8,opt1_fase)
// init_3darray log_parSf(1,2,1,2,1,nbloq_selflo,optSf_fase)

 init_bounded_vector log_A50flom(1,nbloq_selflo,0.67,1.94,optSf_fase);
 init_vector log_sdA50flom(1,nbloq_selflo,optSf_fase);

 init_bounded_vector log_A50floh(1,nbloq_selflo,0.67,1.94,optSf_fase);  
 init_vector log_sdA50floh(1,nbloq_selflo,optSf_fase);

 init_bounded_vector log_A50cruh(1,nbloq_selcru,0.67,1.94,optSc_fase);  
 init_vector log_sdA50cruh(1,nbloq_selcru,optSc_fase);

 init_bounded_vector log_A50crum(1,nbloq_selcru,0.67,1.94,optSc_fase);  
 init_vector log_sdA50crum(1,nbloq_selcru,optSc_fase);


// parametros reclutamientos y mortalidades)
 init_number log_Ro(1)
 init_bounded_number log_pRm(-2.3,-0.1,opt_Rm) // prop de machos en el reclutamiento
 init_bounded_dev_vector dev_log_Ro(1,ntime,-10,10,opt_devRt)
 init_bounded_vector dev_log_Nom(1,nedades,-10,10,opt_devNo)
 init_bounded_vector dev_log_Noh(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_Fm(1,ntime,-20,-0.2,opt_F) // log  mortalidad por pesca por flota
 init_bounded_vector log_Fh(1,ntime,-20,-0.2,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qflo(1,nbloq_qflo,opt_qf)
 init_vector log_qcru(1,nbloq_qcru,opt_qc)

// Crecimimiento
 init_number log_Lom(opt_Lo)
 init_number log_cv_edadm(opt_cva)

 init_number log_Loh(opt_Lo)
 init_number log_cv_edadh(opt_cva)

 init_number log_Mh(opt_M)
 init_number log_Mm(opt_M)

// Fpbr
 init_vector log_Fref(1,npbr,opt_Fpbr);

//---------------------------------------------------------------------------------

//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)

 vector Brec(1,ntime)
// vector pred_CPUE(1,ntime);
// vector pred_Bcru(1,ntime);

// vector pred_Desemb(1,ntime);
 vector likeval(1,20);
 vector Neqm(1,nedades);
 vector Neqh(1,nedades);

 vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_anos(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector mu_edadm(1,nedades)
 vector mu_edadh(1,nedades)
 vector sigma_edadm(1,nedades)
 vector sigma_edadh(1,nedades)
 vector BDo(1,ntime);
 vector No(1,nedades)
 vector prior(1,7)
 vector prop_hpred(1,ntime);

 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Bcru(1,ntime);
 vector prop_h(1,ntime);

 vector Lobs(1,ntime);
 vector Lpred(1,ntime);

 matrix cv_index(1,4,1,ntime)
 matrix nm_sex(1,4,1,ntime)

 matrix Sel_crum(1,ntime,1,nedades);
 matrix Sel_cruh(1,ntime,1,nedades);

 matrix S1(1,nbloq_selflo,1,nedades)
 matrix S2(1,nbloq_selflo,1,nedades)
 matrix S3(1,nbloq_selcru,1,nedades)
 matrix S4(1,nbloq_selcru,1,nedades)


 matrix Sel_m(1,ntime,1,nedades)
 matrix Sel_h(1,ntime,1,nedades)
 matrix Fm(1,ntime,1,nedades)
 matrix Fh(1,ntime,1,nedades)
 matrix Zm(1,ntime,1,nedades)
 matrix Zh(1,ntime,1,nedades)
 matrix Sm(1,ntime,1,nedades)
 matrix Sh(1,ntime,1,nedades)


 matrix Nm(1,ntime,1,nedades)
 matrix Nh(1,ntime,1,nedades)

 matrix NM(1,ntime,1,nedades)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo_m(1,ntime,1,ntallas)
 matrix NVflo_h(1,ntime,1,ntallas)
 matrix NVcru_m(1,ntime,1,ntallas)
 matrix NVcru_h(1,ntime,1,ntallas)

 matrix pred_Ctotm(1,ntime,1,ntallas)
 matrix pred_Ctot_am(1,ntime,1,nedades)
 matrix pred_Ctoth(1,ntime,1,ntallas)
 matrix pred_Ctot_ah(1,ntime,1,nedades)

 matrix pobs_m(1,ntime,1,ntallas)
 matrix ppred_m(1,ntime,1,ntallas)
 matrix pobs_h(1,ntime,1,ntallas)
 matrix ppred_h(1,ntime,1,ntallas)


 matrix pobsc_m(1,ntime,1,ntallas)
 matrix ppredc_m(1,ntime,1,ntallas)
 matrix pobsc_h(1,ntime,1,ntallas)
 matrix ppredc_h(1,ntime,1,ntallas)

 matrix Prob_talla_m(1,nedades,1,ntallas)
 matrix Prob_talla_h(1,nedades,1,ntallas)

 matrix P1(1,nedades,1,ntallas)
 matrix P2(1,nedades,1,ntallas)
 matrix P3(1,nedades,1,ntallas)
 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,nedades)

 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
 number suma6
 number suma7

 number penalty

 number So
 number alfa
 number beta

 number Linfm
 number km
 number Linfh
 number kh

 number Mm
 number Mh

 number BDp
 number Npplus
 number Bph
 number Bpm
 number nm1;
 number cuenta1;
 number nm2;
 number cuenta2;
 number nm3;
 number cuenta3;
 number nm4;
 number cuenta4;

 vector Nph(1,nedades)
 vector Zpbrh(1,nedades)
 vector Fpbrh(1,nedades)
 vector Sph(1,nedades)
 vector Npm(1,nedades)
 vector Zpbrm(1,nedades)
 vector Fpbrm(1,nedades)
 vector Spm(1,nedades)

 vector CTP(1,ntallas)
 vector NMDp(1,ntallas);
 matrix YTP(1,ntime_sim,1,npbr)
  matrix BTp(1,ntime_sim,1,npbr)

 number BD_lp
 vector ratio_pbr(1,npbr)

 vector Nvp(1,nedades);
 number Nvplus;
 vector SDvp(1,ntime_sim);
  
 sdreport_vector CBA(1,npbr) // 

 objective_function_value f
 sdreport_vector pred_CPUE(1,ntime) // 
 sdreport_vector pred_Bcru(1,ntime) // 
 sdreport_vector pred_Desemb(1,ntime) // 
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector BV(1,ntime) 
 sdreport_vector RPR(1,ntime) // 
 
 sdreport_number SSBo
 sdreport_vector RPRp(1,npbr) // RPR proyectado en la simulacion
 sdreport_vector Restim(1,ntime)//Reclutas hembras
 sdreport_vector RPRlp(1,ntime)//
 sdreport_matrix SSBp(1,ntime_sim,1,npbr)// Biomasa desovante proyectada 
 vector Lmf_obs(1,ntime)
 sdreport_vector Lmf_pred(1,ntime)
 vector Lhf_obs(1,ntime)
 sdreport_vector Lhf_pred(1,ntime)
 vector Lmc_obs(1,ntime)
 sdreport_vector Lmc_pred(1,ntime)
 vector Lhc_obs(1,ntime)
 sdreport_vector Lhc_pred(1,ntime)
 sdreport_vector Frpr(1,ntime);


PRELIMINARY_CALCS_SECTION

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


RUNTIME_SECTION
 // maximum_function_evaluations 500,2000,5000
 //convergence_criteria  1e-2,1e-5,1e-5

 // tomado de AMAK (Ianelli)
//  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-5,1e-5

  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-4
  maximum_function_evaluations 100,100,200,300,3500

PROCEDURE_SECTION
// se listan las funciones que contienen los calculos
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



FUNCTION Eval_prob_talla_edad


 int i, j;

// genero una clave edad-talla para otros calculos. Se modela desde L(1)
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

//----------------------------------------------------------------------
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
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
//----------------------------------------------------------------------

FUNCTION Eval_selectividad
 int i,j;


// FLOTA

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

// FLOTA
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


FUNCTION Eval_mortalidades

 Mm=exp(log_Mm);
 Mh=exp(log_Mh);
 
 Fm=elem_prod(Sel_m,outer_prod(mfexp(log_Fm),Unos_edad));
 Fh=elem_prod(Sel_h,outer_prod(mfexp(log_Fh),Unos_edad));

 Zm=Fm+Mm;
 Zh=Fh+Mh;

 Sm=mfexp(-1.0*Zm);
 Sh=mfexp(-1.0*Zh);


FUNCTION Eval_abundancia
 int i, j;


 // Biomasa desovante virgen de largo plazo
 No(1)=exp(log_Ro); // hembras
 for (int j=2;j<=nedades;j++)
 {
   No(j)=No(j-1)*exp(-1.*Mh);
 }
//   No(nedades)+=No(nedades)*exp(-1.*Mh);
   No(nedades)=No(nedades)*exp(-1.*Mh)/(1-exp(-1.*Mh));
  SSBo=sum(elem_prod(No*exp(-dt(1)*Mh)*Prob_talla_h,elem_prod(msex,Wmed(2))));

// Stock-recluta

 alfa=4*h*exp(log_Ro)/(5*h-1);//
 beta=(1-h)*SSBo/(5*h-1);//  

// genero una estructura inicial en equilibrio para el primer año

 Neqh(1)=mfexp(log_Ro);//hembras
 for (j=2;j<=nedades;j++)
 {
   Neqh(j)=Neqh(j-1)*exp(-Zh(1,j-1));
 }
//   Neqh(nedades)+=Neqh(nedades)*exp(-1.*Zh(1,nedades)); // MODIFICAR POR LA OTRA FORMA
   Neqh(nedades)=Neqh(nedades)*exp(-1.*Zh(1,nedades))/(1-exp(-1.*Zh(1,nedades)));



 Neqm(1)=mfexp(log_Ro) * (exp(log_pRm)/(1-exp(log_pRm)));//machos
 for (j=2;j<=nedades;j++)
 {
   Neqm(j)=Neqm(j-1)*exp(-Zm(1,j-1));
 }
   //Neqm(nedades)+=Neqm(nedades)*exp(-1.*Zm(1,nedades));//
   Neqm(nedades)=Neqm(nedades)*exp(-1.*Zm(1,nedades))/(1-exp(-1.*Zm(1,nedades)));



// Abundancia inicial

 Nm(1)=elem_prod(Neqm,exp(dev_log_Nom));
 Nh(1)=elem_prod(Neqh,exp(dev_log_Noh));
 BD(1)=sum(elem_prod(elem_prod(Nh(1),exp(-dt(1)*Zh(1)))*Prob_talla_h,elem_prod(msex,Wmed(2))));

 Rpred(1)=mfexp(log_Ro);//


// se estima la sobrevivencia por edad(a+1) y año(t+1)
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

FUNCTION Eval_deinteres

// Rutina para calcular RPR
 Nv=Nh;// solo para empezar los calculos

// se estima la sobrevivencia por edad(a+1) y año(t+1)
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*Mh);
     Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*Mh);// grupo plus
 }


 NDv=elem_prod((Nv*exp(-dt(1)*Mh))*Prob_talla_h,outer_prod(Unos_anos,msex));
 BDo=NDv*Wmed(2);
 RPR=elem_div(BD,BDo);


 RPRlp=BD/SSBo;


FUNCTION Eval_biomasas
 
 NMD=elem_prod(Nh,mfexp(-dt(1)*Zh))*Prob_talla_h;
 NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));
 
 NVflo_m=elem_prod(elem_prod(Nm,mfexp(-dt(2)*(Zm))),Sel_m)*Prob_talla_m;
 NVflo_h=elem_prod(elem_prod(Nh,mfexp(-dt(2)*(Zh))),Sel_h)*Prob_talla_h;

 NVcru_m=elem_prod(elem_prod(Nm,mfexp(-dt(3)*(Zm))),Sel_crum)*Prob_talla_m;
 NVcru_h=elem_prod(elem_prod(Nh,mfexp(-dt(3)*(Zh))),Sel_cruh)*Prob_talla_h;

// vectores de biomasas derivadas
 BD=NMD*Wmed(2);
 BMflo=NVflo_m*Wmed(1)+NVflo_h*Wmed(2);
 BMcru=NVcru_m*Wmed(1)+NVcru_h*Wmed(2);
 BV=BMflo;

 BT=(Nm*Prob_talla_m)*Wmed(1)+(Nh*Prob_talla_h)*Wmed(2);  


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y año
 pred_Ctot_am=elem_prod(elem_div(Fm,Zm),elem_prod(1.-Sm,Nm));
 pred_Ctotm=pred_Ctot_am*Prob_talla_m;

 pred_Ctot_ah=elem_prod(elem_div(Fh,Zh),elem_prod(1.-Sh,Nh));
 pred_Ctoth=pred_Ctot_ah*Prob_talla_h;

// Proporción total anual de hembras en las capturas
 prop_hpred = elem_div(rowsum(pred_Ctoth),rowsum(pred_Ctoth+pred_Ctotm+1e-10));

// vectores de desembarques predichos por año
 pred_Desemb=pred_Ctotm*Wmed(1)+pred_Ctoth*Wmed(2);


// PROPORCIONES  matrices de proporcion de capturas por talla y año
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


FUNCTION Eval_indices
 

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

FUNCTION Eval_PBR


 for (int i=1;i<=npbr;i++){

 Fpbrh=Sel_h(ntime)*exp(log_Fref(i));
 Zpbrh=Fpbrh+Mh;

 Neqh(1)=mfexp(log_Ro);//hembras
 for (int j=2;j<=nedades;j++)
 {
   Neqh(j)=Neqh(j-1)*exp(-Zpbrh(j-1));
 }
//   Neqh(nedades)+=Neqh(nedades)*exp(-1.*Zpbr(nedades)); // MODIFICAR POR LA OTRA FORMA
   Neqh(nedades)=Neqh(nedades)*exp(-1.*Zpbrh(nedades))/(1-exp(-1.*Zpbrh(nedades))); // MODIFICAR POR LA OTRA FORMA


 BD_lp=sum(elem_prod(elem_prod(Neqh,exp(-dt(1)*Zpbrh))*Prob_talla_h,elem_prod(msex,Wmed(2))));
 ratio_pbr(i)=BD_lp/SSBo;
 }

 Frpr=mfexp(log_Fh)/mfexp(log_Fref(3));
 
FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
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


FUNCTION Eval_funcion_objetivo

 suma4=0; suma5=0; suma6=0; suma7=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//Crucero
 likeval(3)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 likeval(4)=0.5*suma3;// prop p_hembras
 likeval(5)=-1.*sum(nm_sex(1)*elem_prod(pobs_m,log(ppred_m)));
 likeval(6)=-1.*sum(nm_sex(2)*elem_prod(pobs_h,log(ppred_h)));
 likeval(7)=-1.*sum(nm_sex(3)*elem_prod(pobsc_m,log(ppredc_m)));
 likeval(8)=-1.*sum(nm_sex(4)*elem_prod(pobsc_h,log(ppredc_h)));

// lognormal Ninicial y Reclutas
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


//---------------------------------------------------------------------


FUNCTION Eval_CTP

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



REPORT_SECTION

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


// ESTIMA nm y CV

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


FUNCTION Eval_mcmc
  if(reporte_mcmc == 0)
  mcmc_report<<"Bcru_last CTP1 CTP2 CTP3 CTP4 BDp1/BDlast BDp2/BDlast BDp3/BDlast BDp4/BDlast "<<endl;
  mcmc_report<<pred_Bcru(ntime)<<" "<<YTP(2,1)<<" "<<YTP(2,2)<<" "<<YTP(2,3)<<" "<<YTP(2,4)<<
     " "<<SSBp(ntime_sim,1)/BD(ntime)<<" "<<SSBp(ntime_sim,2)/BD(ntime)<<" "<<SSBp(ntime_sim,3)/BD(ntime)<<
     " "<<SSBp(ntime_sim,4)/BD(ntime)<<endl;

  reporte_mcmc++;



FINAL_SECTION

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


