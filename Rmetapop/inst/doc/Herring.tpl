GLOBALS_SECTION
 #include <admodel.h>
 ofstream OutReport;
 
// ======================================================================================================================

DATA_SECTION

 init_int Npop;                                                       // Number of populations
 init_int MaxAge;                                                     // Maximum age
 int InitYr;                                                          // Burnin 
 !! InitYr = -40;
 init_int Nyear;                                                      // Number of years
 init_int Nproj;
 init_int Nfleet;                                                     // Number of fleets
 init_int SRType;                                                     // Type of sr-relationship
 init_int PhaseDiffusionA;
 init_int PhaseDiffusionJ;
 init_int PhaseSigmaR;
 init_int LastRecDev;                                                 // Turn bias correction off
 
 init_4darray CatchInput(1,Npop,1,Nfleet,1,Nyear,-2,1);               // Catches
 init_3darray Index(1,Npop,1,Nyear,-1,2);                             // Index     
 
 init_3darray SpawnAge(1,Npop,1,Nyear,-3,MaxAge);                     // Spawning index data
 init_4darray FishAgeData(1,Npop,1,Nfleet,1,Nyear,-4,MaxAge);         // Spawning index data

 init_number M;                                                       // Natural mortality
 !! cout << M << endl;
 init_vector Mat(0,MaxAge);                                           // Maturity-at-age
 vector MatProp(0,MaxAge);                                            // Proportion maturing
 init_vector Wght(0,MaxAge);                                          // Weight-at-age
 init_vector Fec(0,MaxAge);                                           // Fecundity-at-age
 init_matrix S(1,2,0,MaxAge);                                         // Selectivity-at-age
 init_number InitNStart;                                              // Initial N
 
 int TheYear;
 int PhaseAlpha;                                                      // Alpha parameters [hase
 !! PhaseAlpha = 2;
 !! if (SRType==1) PhaseAlpha = -1;

 init_int TestCodeDat;
 !! if (TestCodeDat != 123456)
 !!  { cout << "Code does not match " << TestCodeDat << endl; exit(1); }
 !! cout << "Completed Data File" << endl;
 
// ======================================================================================================================

PARAMETER_SECTION
 init_number Dummy(-1);
 init_vector LogRbar(1,Npop,1);                                          // Virgin recruitment
 init_number SigmaJ(PhaseDiffusionJ);
 init_number SigmaA(PhaseDiffusionA);
 init_vector Alpha(1,Npop,PhaseAlpha);
 init_bounded_number SigmaR(0,5,PhaseSigmaR);                                                    // SigmaR    
 init_matrix Rec_dev(1,Npop,1,LastRecDev,1);                           // Recruitment deviations

 vector Beta(1,Npop);

 4darray N(1,Npop,InitYr,Nyear+Nproj+1,1,2,0,MaxAge);                  // N-matrix
 4darray N_check(1,Npop,InitYr,Nyear+Nproj+1,1,2,0,MaxAge);            // N_check-matrix
 3darray N1(1,Npop,1,2,0,MaxAge);                                      // After mortality
 3darray N2(1,Npop,1,2,0,MaxAge);                                      // After diffusion
 3darray N3(1,Npop,1,2,0,MaxAge);                                      // After 2nd fishery
 matrix SSB(1,Npop,InitYr,Nyear+Nproj+1);
 matrix SSBP(1,Npop,InitYr,Nyear+Nproj+1);                             // SSB (before fishery)
 matrix Eggs(1,Npop,InitYr,Nyear+Nproj+1);
 3darray F(1,Npop,1,2,InitYr,Nyear+Nproj+1);
 3darray SSBAGE(1,Npop,InitYr,Nyear+Nproj+1,0,MaxAge);                 // Spawning age data
 4darray CAgePred(1,Nfleet,1,Npop,InitYr,Nyear+Nproj,0,MaxAge);        // Predicted catch-at-age
 
 matrix Dist(1,Npop,1,Npop);                                           // Distance between populations
 matrix DiffuseA(1,Npop,1,Npop);                                       // Diffussion (adults)
 matrix DiffuseJ(1,Npop,1,Npop);                                       // Diffussion (juveniles)

 vector LikeIndex(1,Npop);                                             // Index Likelihood
 vector LikeSSBAge(1,Npop);                                            // SSB age data likelihood
 matrix LikeFishAge(1,2,1,Npop);
 number Rec_Penal;                                                     // Recruitment penalty
 number Catch_Penal;
 
 number Fproj;                                                         // Projected F
 matrix Quota(1,Npop,Nyear+1,Nyear+Nproj);                             // Quotas
 

 sdreport_matrix SSBO(1,Npop,InitYr,Nyear+Nproj+1);
 
 
 objective_function_value Obj;

 !! cout << "Completed Parameter Specification" << endl;

// ======================================================================================================================

PRELIMINARY_CALCS_SECTION
 int Iage,Ipop,Jpop,Iprob1,Iprob2,Iprob3;
 
 // Proportion maturing
 MatProp(0) = Mat(0);
 for (Iage=1;Iage<=MaxAge;Iage++)
  if (Mat(Iage-1) >= 1) MatProp(Iage) = 1; else MatProp(Iage) = (Mat(Iage)-Mat(Iage-1))/(1-Mat(Iage-1));
 
 // Distance between populations
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Jpop=1;Jpop<=Npop;Jpop++)
   {
    Iprob1 = abs(Ipop-Jpop);
    Iprob2 = abs(Ipop-(Jpop-Npop));
    Iprob3 = abs((Ipop-Npop)-Jpop);
    Dist(Ipop,Jpop) = min(min(Iprob1,Iprob2),Iprob3);
   }
 

 cout << "Completed Preliminary Calcs Section" << endl;

// ======================================================================================================================
PROCEDURE_SECTION
 int Ipop;
 
 // Specify basic parameters
 SetParameters();
 
 // Initial model state
 InitialState();
 
 // Now project
 Project();
 
 // Likelihood function
 LikeLihood();
 
 Obj = Dummy*Dummy;
 for (Ipop=1;Ipop<=Npop;Ipop++) Obj += LikeIndex(Ipop);
 for (Ipop=1;Ipop<=Npop;Ipop++) Obj += LikeSSBAge(Ipop);
 for (Ipop=1;Ipop<=Npop;Ipop++) Obj += (LikeFishAge(1,Ipop)+LikeFishAge(2,Ipop));
 Obj += (Rec_Penal + Catch_Penal);
 
 //cout << LikeIndex << endl;
 //cout << LikeSSBAge << endl;
 //cout << LikeFishAge << endl;
 //cout << Rec_Penal << " " << Catch_Penal << endl;
 
 // Store SSB
 SSBO = SSB;
 
 cout << Obj << endl;
 
// ======================================================================================================================

FUNCTION SetParameters
 int Ipop,Jpop;
 dvariable TotalJ, TotalA;
 
 for (Ipop=1;Ipop<=Npop;Ipop++)
  {
   TotalJ = 0; TotalA = 0;
   for (Jpop=1;Jpop<=Npop;Jpop++)
    {
     DiffuseJ(Ipop,Jpop) = exp(-1*square(Dist(Ipop,Jpop))/square(SigmaJ));
     TotalJ += DiffuseJ(Ipop,Jpop);
     DiffuseA(Ipop,Jpop) = exp(-1*square(Dist(Ipop,Jpop))/square(SigmaA));
     TotalA +=DiffuseA(Ipop,Jpop);
    }
   for (Jpop=1;Jpop<=Npop;Jpop++)
    {
     DiffuseJ(Ipop,Jpop) /= TotalJ;
     DiffuseA(Ipop,Jpop) /= TotalA;
    }
  }  
 
// ======================================================================================================================

FUNCTION InitialState
 int Iage,Ipop,Jpop,Imat,Iyear;
 dvariable Pmat,Temp;
 dvar_vector NextN(1,Npop);
 dvar_matrix MatCome(1,Npop,0,MaxAge);

 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Iage=0;Iage<=MaxAge;Iage++)
   for (Imat=1;Imat<=2;Imat++)
    N(Ipop,InitYr,Imat,Iage) = InitNStart;

 for (Iyear=InitYr;Iyear<=0;Iyear++)
  {
   // Remove natural mortality
   for (Ipop=1;Ipop<=Npop;Ipop++)
    for (Iage=0;Iage<=MaxAge;Iage++)
     for (Imat=1;Imat<=2;Imat++)
      N1(Ipop,Imat,Iage) = N(Ipop,Iyear,Imat,Iage)*exp(-M);
      
   // Difussion (1+)
   for (Iage=1;Iage<=MaxAge;Iage++)
    for (Imat=1;Imat<=2;Imat++)
     {
      NextN.initialize();
      for (Ipop=1;Ipop<=Npop;Ipop++)
       for (Jpop=1;Jpop<=Npop;Jpop++)
        NextN(Jpop) += DiffuseA(Ipop,Jpop)*N1(Ipop,Imat,Iage);
      for (Ipop=1;Ipop<=Npop;Ipop++)
       N2(Ipop,Imat,Iage) = NextN(Ipop);      
     }
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Imat=1;Imat<=2;Imat++)
      N2(Ipop,Imat,0) = N1(Ipop,Imat,0);
     
    // Maturation
    MatCome.initialize();
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Iage=1;Iage<=MaxAge;Iage++)
      {
       // Numbers maturing
       Pmat = MatProp(Iage) * N2(Ipop,1,Iage);
       for (Jpop=1;Jpop<=Npop;Jpop++)
        MatCome(Jpop,Iage) += Pmat*DiffuseJ(Ipop,Jpop);
      }
    
    // Account for maturation
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Iage=1;Iage<=MaxAge;Iage++)
      {
       N3(Ipop,2,Iage) = N2(Ipop,2,Iage) + MatCome(Ipop,Iage);
       N3(Ipop,1,Iage) = N2(Ipop,1,Iage)*(1-MatProp(Iage));
      }
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Imat=1;Imat<=2;Imat++)
      N3(Ipop,Imat,0) = N2(Ipop,Imat,0);
      
    // Update dynamics
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Imat=1;Imat<=2;Imat++)
      {
       for (Iage=0;Iage<MaxAge;Iage++) N(Ipop,Iyear+1,Imat,Iage+1) = N3(Ipop,Imat,Iage);
       N(Ipop,Iyear+1,Imat,MaxAge) += N3(Ipop,Imat,MaxAge);
      } 
       
    // SSB
    for (Ipop=1;Ipop<=Npop;Ipop++)
     {
      SSB(Ipop,Iyear) = 0;
      for (Iage=1;Iage<=MaxAge;Iage++) SSB(Ipop,Iyear) += N3(Ipop,2,Iage)*Wght(Iage);
      SSBP(Ipop,Iyear) = SSB(Ipop,Iyear);
     }
     
    // Recruitment (all immature)
    for (Ipop=1;Ipop<=Npop;Ipop++)
     {
      N(Ipop,Iyear+1,1,0) = mfexp(LogRbar(Ipop));
      N(Ipop,Iyear+1,2,0) = 0;
     } 
  }
 for (Ipop=1;Ipop<=Npop;Ipop++) 
  { 
   if (SRType==2)                                        // Ricker
    {
     Temp = N(Ipop,0,1,0)/(Alpha(Ipop)*SSB(Ipop,0));
     Temp = -1.0*log(Temp); 
     Beta(Ipop) = Temp/SSB(Ipop,0);
    } 
   if (SRType==3)                                        // Beverton-Holt
    {
     Temp = (Alpha(Ipop)*SSB(Ipop,0))/N(Ipop,0,1,0) - 1.0;
     Beta(Ipop) = Temp/SSB(Ipop,0);
    }
  } 
  
// ======================================================================================================================

FUNCTION ProjectOneYear
 int Iage,Ipop,Jpop,Imat,Iyear,tune_F,F_tune,max_harvest_rate;
 dvariable Pmat,vbio,temp,join1,temp1,Z_adjuster,Z_adjuster2,RMult,ExpRec;
 dvar_vector NextN(1,Npop);
 dvar_matrix MatCome(1,Npop,0,MaxAge);
 dvar_matrix Z_rate(1,2,0,MaxAge),Z_rate2(1,2,0,MaxAge);

 // Set the year
 Iyear = TheYear;

 // Set F_tune and the maximum harvest rate
 F_tune = 5;
 max_harvest_rate = 3;

 // Solve for catch
 for (Ipop=1;Ipop<=Npop;Ipop++)
  {
   if (Iyear <= Nyear)
    if (CatchInput(Ipop,1,Iyear,1) > 0)
     { 
      // Initial guess
      vbio = 0;
      for (Iage=0;Iage<=MaxAge;Iage++)
       for (Imat=1;Imat<=2;Imat++)
        vbio += N(Ipop,Iyear,Imat,Iage)*S(Imat,Iage)*Wght(Iage);     
      temp = CatchInput(Ipop,1,Iyear,1)/(vbio + CatchInput(Ipop,1,Iyear,1));    
      join1=1.0/(1.0+mfexp(30.*(temp-0.95)));
      temp1=join1*temp + (1.0-join1)*0.95;
      F(Ipop,1,Iyear) = -log(1.-temp1);
      
      // Tune the Fs
      for (tune_F=1;tune_F<=F_tune;tune_F++)
       {
        // Set up Z
        for (Iage=0;Iage<=MaxAge;Iage++)
         for (Imat=1;Imat<=2;Imat++)
          {
           Z_rate(Imat,Iage) = M + S(Imat,Iage)*F(Ipop,1,Iyear);
           Z_rate2(Imat,Iage) = (1-mfexp(-Z_rate(Imat,Iage)))/Z_rate(Imat,Iage);  
          }         
        
        // Now tune
        if (tune_F < F_tune)
         {
          Z_adjuster2 = 0;
          for (Iage=0;Iage<=MaxAge;Iage++)
           for (Imat=1;Imat<=2;Imat++)
            Z_adjuster2 += F(Ipop,1,Iyear)*N(Ipop,Iyear,Imat,Iage)*S(Imat,Iage)*Wght(Iage)*Z_rate2(Imat,Iage);
          Z_adjuster = CatchInput(Ipop,1,Iyear,1)/(Z_adjuster2+0.0001);
       
          // Adjust total Z
          for (Iage=0;Iage<=MaxAge;Iage++)
           for (Imat=1;Imat<=2;Imat++)
            {
             Z_rate(Imat,Iage)  = M + Z_adjuster*(Z_rate(Imat,Iage)-M);
             Z_rate2(Imat,Iage) = (1-mfexp(-Z_rate(Imat,Iage)))/Z_rate(Imat,Iage);  
            }
        
          // Adjust total exploitable biomass
          Z_adjuster2 = 0;
          for (Iage=0;Iage<=MaxAge;Iage++)
           {
            CAgePred(1,Ipop,Iyear,Iage) = 0;
            for (Imat=1;Imat<=2;Imat++)
             {
              Z_adjuster2 += N(Ipop,Iyear,Imat,Iage)*S(Imat,Iage)*Wght(Iage)*Z_rate2(Imat,Iage);
              CAgePred(1,Ipop,Iyear,Iage) += F(Ipop,1,Iyear)*N(Ipop,Iyear,Imat,Iage)*S(Imat,Iage)*Z_rate2(Imat,Iage);
             } 
           }  
          temp = CatchInput(Ipop,1,Iyear,1)/(Z_adjuster2 + 0.00001);    
          join1=1.0/(1.0+mfexp(30.*(temp-0.95*max_harvest_rate)));
          F(Ipop,1,Iyear) = join1*temp + (1.0-join1)*max_harvest_rate;
         }
       }
      //cout << "A" << Ipop << " " << Iyear << " " << F(Ipop,1,Iyear) << " " << CatchInput(Ipop,1,Iyear,1) << " " << F(Ipop,1,Iyear)*Z_adjuster2 << " " << Z_adjuster2 << endl;
      Catch_Penal += square(CatchInput(Ipop,1,Iyear,1)-F(Ipop,1,Iyear)*Z_adjuster2);
     }
    else
     {
      F(Ipop,1,Iyear) = 0;
      for (Iage=0;Iage<=MaxAge;Iage++)
       for (Imat=1;Imat<=2;Imat++)
        Z_rate(Imat,Iage)  = M;
     }
   else
     {
      F(Ipop,1,Iyear) = 0;
      for (Iage=0;Iage<=MaxAge;Iage++)
       for (Imat=1;Imat<=2;Imat++)
        Z_rate(Imat,Iage)  = M;
     }
    
   for (Iage=0;Iage<=MaxAge;Iage++)
    for (Imat=1;Imat<=2;Imat++)
     N1(Ipop,Imat,Iage) = N(Ipop,Iyear,Imat,Iage)*exp(-Z_rate(Imat,Iage));
  }
      
 // Difussion (1+)
 for (Iage=1;Iage<=MaxAge;Iage++)
  for (Imat=1;Imat<=2;Imat++)
   {
    NextN.initialize();
    for (Ipop=1;Ipop<=Npop;Ipop++)
     for (Jpop=1;Jpop<=Npop;Jpop++)
      NextN(Jpop) += DiffuseA(Ipop,Jpop)*N1(Ipop,Imat,Iage);
    for (Ipop=1;Ipop<=Npop;Ipop++)
     N2(Ipop,Imat,Iage) = NextN(Ipop);      
   }
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Imat=1;Imat<=2;Imat++)
    N2(Ipop,Imat,0) = N1(Ipop,Imat,0);
     
  // Maturation
  MatCome.initialize();
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Iage=1;Iage<=MaxAge;Iage++)
    {
     // Numbers maturing
     Pmat = MatProp(Iage) * N2(Ipop,1,Iage);
     for (Jpop=1;Jpop<=Npop;Jpop++)
      MatCome(Jpop,Iage) += Pmat*DiffuseJ(Ipop,Jpop);
    }
    
  // Account for maturation
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Iage=1;Iage<=MaxAge;Iage++)
    {
     N3(Ipop,2,Iage) = N2(Ipop,2,Iage) + MatCome(Ipop,Iage);
     N3(Ipop,1,Iage) = N2(Ipop,1,Iage)*(1-MatProp(Iage));
    }
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Imat=1;Imat<=2;Imat++)
    N3(Ipop,Imat,0) = N2(Ipop,Imat,0);

  // Compute pre-fishery SSB
  for (Ipop=1;Ipop<=Npop;Ipop++)
   {
    SSBP(Ipop,Iyear) = 0;
    for (Iage=1;Iage<=MaxAge;Iage++) 
     SSBP(Ipop,Iyear) += N3(Ipop,2,Iage)*Wght(Iage);
   }

  // Remove 2nd fishery catch
  for (Ipop=1;Ipop<=Npop;Ipop++)
   {
    if (Iyear <= Nyear)
     {
      if (CatchInput(Ipop,2,Iyear,1) > 0)
       { 
        // Initial guess
        vbio = 0;
        for (Iage=0;Iage<=MaxAge;Iage++) vbio += N3(Ipop,2,Iage)*Wght(Iage);     
        temp = CatchInput(Ipop,2,Iyear,1)/(vbio + CatchInput(Ipop,2,Iyear,1));    
        join1=1.0/(1.0+mfexp(30.*(temp-0.95)));
        temp1=join1*temp + (1.0-join1)*0.95;
        F(Ipop,2,Iyear) = -log(1.-temp1);
           
        // Tune the Fs
        for (tune_F=1;tune_F<=F_tune;tune_F++)
         {
          // Set up Z
          for (Iage=0;Iage<=MaxAge;Iage++)
           {
            Z_rate(2,Iage) = F(Ipop,2,Iyear);
            Z_rate2(2,Iage) = (1-mfexp(-Z_rate(2,Iage)))/Z_rate(2,Iage);  
           }         
              
          // Now tune
          if (tune_F < F_tune)
           {
            Z_adjuster2 = 0;
            for (Iage=0;Iage<=MaxAge;Iage++) Z_adjuster2 += F(Ipop,2,Iyear)*N3(Ipop,2,Iage)*Wght(Iage)*Z_rate2(2,Iage);
            Z_adjuster = CatchInput(Ipop,2,Iyear,1)/(Z_adjuster2+0.0001);
             
            // Adjust total Z
            for (Iage=0;Iage<=MaxAge;Iage++)
             {
              Z_rate(2,Iage)  = Z_adjuster*Z_rate(2,Iage);
              Z_rate2(2,Iage) = (1-mfexp(-Z_rate(2,Iage)))/Z_rate(2,Iage);  
             }
             
            // Adjust total exploitable biomass
            Z_adjuster2 = 0;
            for (Iage=0;Iage<=MaxAge;Iage++) 
             {
              Z_adjuster2 += N3(Ipop,2,Iage)*Wght(Iage)*Z_rate2(2,Iage);
              CAgePred(2,Ipop,Iyear,Iage) = N3(Ipop,2,Iage)*Z_rate2(2,Iage);
             } 
            temp = CatchInput(Ipop,2,Iyear,1)/(Z_adjuster2 + 0.00001);    
            join1=1.0/(1.0+mfexp(30.*(temp-0.95*max_harvest_rate)));
            F(Ipop,2,Iyear) = join1*temp + (1.0-join1)*max_harvest_rate;
           }
         }
        //cout << "B" << Ipop << " " << Iyear << " " << F(Ipop,2,Iyear) << " " << CatchInput(Ipop,2,Iyear,1) << " " << F(Ipop,2,Iyear)*Z_adjuster2 << endl;
        Catch_Penal += square(CatchInput(Ipop,2,Iyear,1)-F(Ipop,2,Iyear)*Z_adjuster2);
       }
      else
       F(Ipop,2,Iyear) = 0;
      }    
    else
     {
      F(Ipop,2,Iyear) = Fproj;
      Quota(Ipop,Iyear) = 0;
      for (Iage=0;Iage<=MaxAge;Iage++) 
       {
        Z_rate(2,Iage) = F(Ipop,2,Iyear);
        if (Fproj <= 0)
         Z_rate2(2,Iage) = 0; 
        else 
         Z_rate2(2,Iage) = (1-mfexp(-Z_rate(2,Iage)))/Z_rate(2,Iage);  
        Quota(Ipop,Iyear) += F(Ipop,2,Iyear)*N3(Ipop,2,Iage)*Wght(Iage)*Z_rate2(2,Iage);
       }  
     }
   }
   
  // Remove mature F
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Iage=0;Iage<=MaxAge;Iage++)
    N3(Ipop,2,Iage) *= exp(-F(Ipop,2,Iyear));
    
  // Update dynamics
  for (Ipop=1;Ipop<=Npop;Ipop++)
   for (Imat=1;Imat<=2;Imat++)
    {
     for (Iage=0;Iage<MaxAge;Iage++) N(Ipop,Iyear+1,Imat,Iage+1) = N3(Ipop,Imat,Iage);
     N(Ipop,Iyear+1,Imat,MaxAge) += N3(Ipop,Imat,MaxAge);
     for (Iage=0;Iage<=MaxAge;Iage++) N_check(Ipop,Iyear,Imat,Iage) = N3(Ipop,Imat,Iage);
    } 
       
  // SSB
  for (Ipop=1;Ipop<=Npop;Ipop++)
   {
    SSB(Ipop,Iyear) = 0;
    Eggs(Ipop,Iyear) = 0;
    for (Iage=1;Iage<=MaxAge;Iage++) 
     {
      SSB(Ipop,Iyear) += N3(Ipop,2,Iage)*Wght(Iage);
      Eggs(Ipop,Iyear) += N3(Ipop,2,Iage)*Fec(Iage);
      SSBAGE(Ipop,Iyear,Iage) = N3(Ipop,2,Iage);
     } 
   }
     
  // Recruitment (all immature)
  for (Ipop=1;Ipop<=Npop;Ipop++)
   {
    if (Iyear <= LastRecDev)
     RMult = exp(Rec_dev(Ipop,Iyear)-square(SigmaR)/2.0);
    else
     RMult = 1;
    if (SRType==1) ExpRec = mfexp(LogRbar(Ipop));
    if (SRType==2) ExpRec = Alpha(Ipop)*Eggs(Ipop,Iyear)*mfexp(-Beta(Ipop)*Eggs(Ipop,Iyear));
    if (SRType==3) ExpRec = Alpha(Ipop)*Eggs(Ipop,Iyear)/(1.0+Beta(Ipop)*Eggs(Ipop,Iyear));
    N(Ipop,Iyear+1,1,0) = ExpRec*RMult;
     
    N(Ipop,Iyear+1,2,0) = 0;
   } 


// ======================================================================================================================

FUNCTION Project

 F.initialize();

 // Now project with catches
 Catch_Penal = 0;
 for (TheYear=1;TheYear<=Nyear;TheYear++)
  {
   ProjectOneYear();
  }

//=========================================================================================

FUNCTION LikeLihood
 int Ipop,Iyear,Iage,Ifleet;
 dvariable Error,Total1,Total2,Penal;

 // Index likelihood
 LikeIndex.initialize();
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Iyear=1;Iyear<=Nyear;Iyear++)
   if (Index(Ipop,Iyear,1) > 0)
    {
     Error = log(Index(Ipop,Iyear,1)) - log(SSB(Ipop,Iyear));
     //cout << Ipop << " " << Iyear << " " << Index(Ipop,Iyear,1) << " " << SSB(Ipop,Iyear) << endl;
     LikeIndex(Ipop) += 0.5*square(Error) / square(Index(Ipop,Iyear,2));
    }
   
 // SSB age-composition
 LikeSSBAge.initialize();
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Iyear=1;Iyear<=Nyear;Iyear++)
   {
    Total1 = 0;
    for (Iage=1;Iage<=MaxAge;Iage++) Total1 += SpawnAge(Ipop,Iyear,Iage);
    Total2 = 0;
    for (Iage=1;Iage<=MaxAge;Iage++) Total2 += SSBAGE(Ipop,Iyear,Iage);
    for (Iage=1;Iage<=MaxAge;Iage++)
     if (SpawnAge(Ipop,Iyear,Iage) > 0)
      {
       Error = ((SpawnAge(Ipop,Iyear,Iage)/Total1)/(SSBAGE(Ipop,Iyear,Iage)/Total2));
       LikeSSBAge(Ipop) += SpawnAge(Ipop,Iyear,Iage)*log(Error);
      } 
   }
 
 // Catch age-composition
 LikeFishAge.initialize();
 Penal = 0;
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Iyear=1;Iyear<=Nyear;Iyear++)
   for (Ifleet=1;Ifleet<=2;Ifleet++)
    {
     Total1 = 0;
     for (Iage=1;Iage<=MaxAge;Iage++) Total1 += FishAgeData(Ipop,Ifleet,Iyear,Iage);
     Total2 = 0;
     for (Iage=1;Iage<=MaxAge;Iage++) 
      {
       CAgePred(Ifleet,Ipop,Iyear,Iage) = posfun(CAgePred(Ifleet,Ipop,Iyear,Iage),1.0e-10,Penal);
       Total2 += CAgePred(Ifleet,Ipop,Iyear,Iage);
      } 
     for (Iage=1;Iage<=MaxAge;Iage++)
      if (FishAgeData(Ipop,Ifleet,Iyear,Iage) > 0)
       {
        Error = ((FishAgeData(Ipop,Ifleet,Iyear,Iage)/Total1)/(CAgePred(Ifleet,Ipop,Iyear,Iage)/Total2));
        LikeFishAge(Ifleet,Ipop) += FishAgeData(Ipop,Ifleet,Iyear,Iage)*log(Error);
        //cout << " " << Ifleet << " " << Ipop << " " << Iyear << " " << Iage << " " << FishAgeData(Ipop,Ifleet,Iyear,Iage) << " " << CAgePred(Ifleet,Ipop,Iyear,Iage) << " " << Error << " " << log(Error) << endl;
       } 
    }
    
 // Recruitment penalty
 Rec_Penal = 0;
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Iyear=1;Iyear<=LastRecDev;Iyear++)
   Rec_Penal += log(SigmaR) + square(Rec_dev(Ipop,Iyear))/square(SigmaR);
  
  
//=========================================================================================

FUNCTION OutputStuff

 int Ipop,Imat,Iyear,Iage;

 OutReport.close();
 OutReport.open("OutputFile1.Out");

 OutReport<< "Likelihood components" << endl;
 OutReport<< Rec_Penal << endl;
 OutReport<< LikeIndex << " # Log-likelihood for the Index" << endl;
 OutReport<< LikeSSBAge << " # Log-likelihood for the spawning index" << endl;
 OutReport<< LikeFishAge << " # Log-likelihood for the fishery" << endl;

 OutReport<< "#MatProp" << endl; 
 OutReport << MatProp << endl;
 OutReport<< "#Dist" << endl;  
 OutReport << Dist << endl;  
 OutReport<< "#DiffuseJ" << endl;  
 OutReport << DiffuseJ << endl;   
 OutReport<< "#DiffuseA" << endl;   
 OutReport << DiffuseA << endl;   
 OutReport<< "#SigmaA" << endl;   
 OutReport << SigmaA << endl;
 OutReport<< "#SigmaJ" << endl;  
 OutReport << SigmaJ << endl;
 OutReport<< "#mfexpLogRbar" << endl;  
 OutReport << mfexp(LogRbar) << endl;
 OutReport<< "#Alpha" << endl;  
 OutReport << Alpha << endl;
 OutReport<< "#Beta" << endl;  
 OutReport << Beta << endl;
 OutReport<< "#SigmaR" << endl;  
 OutReport << SigmaR << endl;

 OutReport<< "#Quota.Calculation" << endl; 
 Quota.initialize();
 for (Ipop=1;Ipop<=Npop;Ipop++) 
  {
   OutReport << Ipop << " " << SSB(Ipop,Nyear) << " " << SSB(Ipop,0) << " ";
   for (Iyear=Nyear+1;Iyear<=Nyear+Nproj;Iyear++)
    {
     Fproj = 0.000;
     TheYear = Iyear; 
     ProjectOneYear(); 
     if (SSB(Ipop,Iyear) > (0.25*SSB(Ipop,0)*1.25))
      Fproj = 0.2;
     else
     if (SSB(Ipop,Iyear) > (0.25*SSB(Ipop,0)/0.9))
      Fproj =.1;
     else
      Fproj = 0.000;
     TheYear = Iyear; 
     ProjectOneYear(); 
     OutReport << Fproj << " ";
     OutReport << Quota(Ipop,Iyear) << " ";
     OutReport << SSB(Ipop,Iyear) << " ";
    } 
  } 
 OutReport << endl; 


 OutReport<< "#Virgin.conditions" << endl;
 for (Ipop=1;Ipop<=Npop;Ipop++)OutReport<< N(Ipop,0,1,0) << " "; OutReport<< endl;
 for (Ipop=1;Ipop<=Npop;Ipop++)OutReport<< SSB(Ipop,0) << " "; OutReport<< endl;
 for (Ipop=1;Ipop<=Npop;Ipop++)OutReport<< SSBP(Ipop,0) << " "; OutReport<< endl;

 OutReport<< "#For.All.Output.Ages" << endl;
 for (Iyear=InitYr;Iyear<=Nyear;Iyear++)
  {
   OutReport<< Iyear << " ";
   for (Ipop=1;Ipop<=Npop;Ipop++) OutReport<< SSB(Ipop,Iyear) << " " << SSBP(Ipop,Iyear) << " " << N(Ipop,Iyear,1,0) << " ";
   OutReport<< endl;
  }
 
 OutReport<< "#For.All.Output" << endl;
 for (Ipop=1;Ipop<=Npop;Ipop++)
  for (Imat=1;Imat<=2;Imat++)
   for (Iyear=InitYr;Iyear<=Nyear;Iyear++)
    {
     OutReport<< Ipop << " " << Imat << " " << Iyear << " " << SSB(Ipop,Iyear) << " ";
     OutReport<< F(Ipop,1,Iyear) << " " << F(Ipop,2,Iyear) << " ";
     for (Iage=0;Iage<=MaxAge;Iage++)OutReport<< N_check(Ipop,Iyear,Imat,Iage) << " ";
     OutReport<< endl;
    } 

//=========================================================================================

REPORT_SECTION
 OutputStuff();





