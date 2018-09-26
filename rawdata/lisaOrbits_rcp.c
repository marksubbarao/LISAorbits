// ******************************************************************
//                                                                
//MODULE  : lisaOribts.c                                          
//MODIFIED: 16 August 2005                                        
//                                                                
//AUTHORS : Shane L. Larson
//             s.larson@northwestern.edu
//                                                                
//REFERENCES:                                                     
//                                                                
//DESCRIPTION                                                     
//----------------------------------------------------------------
// This program is designed to deduce the standard orbital        
// parameters of a LISA spacecraft orbit from the analytic model  
// of the barycentric coordinate positions for a spacecraft by    
// Rubbo, Cornish and Poujade.  The orbital parameters then will  
// be used in our meteor stream work.                             
//     Rubbo, Cornish & Poujade                                   
//     PRD _69_ 082003 [2004]                                     
//                                                                
// ******************************************************************


//============ INCLUDE LIBRARIES ============

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



//============ CONSTANTS DEFINED ============

#define MSun            1.989E30                   //solar mass in kg  
#define AU              149597900000.0             //m in an AU        
#define lyr             9.46052840488E15           //m in a lyr        
#define pc              3.08567818585E16           //m in a pc         
#define YR              31556925.9747              //seconds per year   

#define PI              3.14159265358979323846     //Why ask pi?       
#define degrees         1.74532925199E-2           //radians per degree

#define G               6.6726E-11                 //Newton G, mks      
#define C               299792458.0                //Speed of light, mks



//============ FUNCTION PROTOTYPES ============

void GetNode(int kASCEND, double time, double delT, double Rearth, double ecc,
     double Torb, double kg, double beta, double *nodeTime, double *nodeAnomaly,
     double *nodeR, double *Omega);

void GetExtremum(int kPERI, double time, double delT, double Rearth, double ecc,
             double Torb, double kg, double beta, double Rvin, double *extrmTime,
             double *extrmAnomaly, double *extrmR);

void GetSpacecraftLocation(double baryctr[],double a, double b, double e, double R);

double vDOT(double v1[], double v2[]);

void Die(char routine[], char errorMsg[]);


//=========================================================
//  MAIN  ROUTINE                                          
//=========================================================

int main(void)
{
   //============= VARIABLE DECLARATIONS =============
   
   unsigned long int nPTS;
   unsigned long int ii, jj;   //looping indexes
   
   double delT, t;
   double kc, kg, alpha, beta1, beta2, beta3;
   double Rearth, ecc, inc, Torb;
   
   double N12[3], N23[3], N13[3];
   double L1, L2, L3;
   double ang1, ang2, ang3;
   
   double s1[3], s2[3], s3[3], s1old[3], s2old[3], s3old[3];
   double th1, th2, th3;
   double R1, R2, R3, R1old, R2old, R3old;
   double drdt1, drdt2, drdt3, drdt1old, drdt2old, drdt3old;
   
   double Tanode1, Tanode2, Tanode3, Tdnode1, Tdnode2, Tdnode3;
   double THanode1, THanode2, THanode3, THdnode1, THdnode2, THdnode3;
   
   double Tapo1, Tapo2, Tapo3, Tperi1, Tperi2, Tperi3;
   double THapo1, THapo2, THapo3, THperi1, THperi2, THperi3;
   
   double Ranode1, Ranode2, Ranode3, Rdnode1, Rdnode2, Rdnode3;
   double Rapo1, Rapo2, Rapo3, Rperi1, Rperi2, Rperi3;
   
   double Omega1, Omega2, Omega3;
   double w1, w2, w3;
   
   double Tmp1;
   
   FILE *ahandle, *s1handle, *s2handle, *s3handle;


   //=============== START ROUTINE HERE ===============

   printf("Welcome to LISA METEOR DELUXE.\n\n");
   
   //==================================================================================
   //---    LISA ORBIT DETERMINAIONS   ------------------------------------------------
   //==================================================================================
   
     
     //---- Time and orbit setup --------------
     
     delT = 5*86400.0;  //spacing of points along orbit in seconds
//    delT = 100.0;  //spacing of points along orbit in seconds
     t = -delT;     //need one point early to set up initial derivatives
     
     nPTS = (unsigned long int)ceil(1.0*YR/delT);
     
     
     kc = 0.0;  //RCP lambda
     kg = 0.0;  //RCP kappa 
     
     
     ahandle = fopen("constellation.dat","w");
     s1handle = fopen("spacecraft1.dat","w");
     s2handle = fopen("spacecraft2.dat","w");
     s3handle = fopen("spacecraft3.dat","w");
     
     fprintf(ahandle,"# idx\ttime(s)\tL1(m)\tL2(m)\tL3(m)\tth1(d)\tth2(d)\tth3(d)\tR1(AU)\tR2(AU)\tR3(AU)\n");
     fprintf(s1handle,"# idx\ttime(s)\tANOMALY(d)\tR1(AU)\ts1x(m)\ts1y(m)\ts1z(m)\n");
     fprintf(s2handle,"# idx\ttime(s)\tANOMALY(d)\tR2(AU)\ts2x(m)\ts2y(m)\ts2z(m)\n");
     fprintf(s3handle,"# idx\ttime(s)\tANOMALY(d)\tR3(AU)\ts3x(m)\ts3y(m)\ts3z(m)\n");
     
     //---- SPACECRAFT CONSTELLATION SETUP -----------------------------
     //-----------------------------------------------------------------
     Rearth = 1.0*AU;
     ecc = (5.0e9)/(2.0*Rearth*sqrt(3.0));
     inc = sqrt(3.0)*ecc;
     
     Torb = sqrt(4.0*PI*PI*Rearth*Rearth*Rearth/(G*MSun));
      
      
     //--- Set up spacecraft at initial timestep -----------
  
     alpha = 2.0*PI*t/Torb + kg;   //RCP alpha
     
     //--- s/c 1 ----
     beta1 = kc;
     GetSpacecraftLocation(s1,alpha,beta1,ecc,Rearth);
     Tmp1 = (alpha - beta1);
     th1 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     //--- s/c 2 ----
     beta2 = 2.0*PI*1.0/3.0 + kc;
     GetSpacecraftLocation(s2,alpha,beta2,ecc,Rearth);
     Tmp1 = (alpha - beta2);
     th2 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     //--- s/c 3 ----
     beta3 = 2.0*PI*2.0/3.0 + kc;
     GetSpacecraftLocation(s3,alpha,beta3,ecc,Rearth);
     Tmp1 = (alpha - beta3);
     th3 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     R1old = sqrt(vDOT(s1,s1));
     R2old = sqrt(vDOT(s2,s2));
     R3old = sqrt(vDOT(s3,s3));
     
     //--- Set up spacecraft at next timestep so derivs can be computed --
     t = 0.0;
     alpha = 2.0*PI*t/Torb + kg;   //RCP alpha
     
     //--- s/c 1 ----
     beta1 = kc;
     GetSpacecraftLocation(s1,alpha,beta1,ecc,Rearth);
     Tmp1 = (alpha - beta1);
     th1 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     //--- s/c 2 ----
     beta2 = 2.0*PI*1.0/3.0 + kc;
     GetSpacecraftLocation(s2,alpha,beta2,ecc,Rearth);
     Tmp1 = (alpha - beta2);
     th2 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     //--- s/c 3 ----
     beta3 = 2.0*PI*2.0/3.0 + kc;
     GetSpacecraftLocation(s3,alpha,beta3,ecc,Rearth);
     Tmp1 = (alpha - beta3);
     th3 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     R1 = sqrt(vDOT(s1,s1));
     R2 = sqrt(vDOT(s2,s2));
     R3 = sqrt(vDOT(s3,s3));
     
     drdt1 = (R1 - R1old)/delT;
     drdt2 = (R2 - R2old)/delT;
     drdt3 = (R3 - R3old)/delT;
     
     
     //-----------------------------------------------
     //  MAIN LOOP OVER LISA ORBIT                    
     //-----------------------------------------------
     
     for (jj = 1; jj < nPTS; jj++)
     {
        t += delT;  //advance timestep
        
        //store point before advancing position
        //used to compute velocity vector      
        for(ii = 0; ii < 3; ii++)
        {
           s1old[ii] = s1[ii];
           s2old[ii] = s2[ii];
           s3old[ii] = s3[ii];
        }
        
        R1old = R1;
        R2old = R2;
        R3old = R3;
        
        drdt1old = drdt1;
        drdt2old = drdt2;
        drdt3old = drdt3;
        

        //get new spacecraft positions & true anomalies ----------
        alpha = 2.0*PI*t/Torb + kg;   //RCP alpha
     
        //--- s/c 1 ----
        beta1 = kc;
        GetSpacecraftLocation(s1,alpha,beta1,ecc,Rearth);
        Tmp1 = (alpha - beta1);
        th1 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
        //--- s/c 2 ----
        beta2 = 2.0*PI*1.0/3.0 + kc;
        GetSpacecraftLocation(s2,alpha,beta2,ecc,Rearth);
        Tmp1 = (alpha - beta2);
        th2 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
        //--- s/c 3 ----
        beta3 = 2.0*PI*2.0/3.0 + kc;
        GetSpacecraftLocation(s3,alpha,beta3,ecc,Rearth);
        Tmp1 = (alpha - beta3);
        th3 = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
        
        R1 = sqrt(vDOT(s1,s1));
        R2 = sqrt(vDOT(s2,s2));
        R3 = sqrt(vDOT(s3,s3));
        
        drdt1 = (R1 - R1old)/delT;
        drdt2 = (R2 - R2old)/delT;
        drdt3 = (R3 - R3old)/delT;
        
        
        //-- construct arm vectors ---
        for (ii = 0; ii < 3; ii++)
        {
           N12[ii] = s2[ii] - s1[ii];
        
           N13[ii] = s3[ii] - s1[ii];

           N23[ii] = s3[ii] - s2[ii];
        }
        
        //normalize arm vectors to unity
        L1 = sqrt(vDOT(N23,N23));
        L2 = sqrt(vDOT(N13,N13));
        L3 = sqrt(vDOT(N12,N12));
		
		Tmp1 = vDOT(N13,N12)/(L2*L3);
		ang1 = acos(Tmp1)*(180.0/PI);
        
		Tmp1 = vDOT(N23,N12)/(L1*L3);
		ang2 = acos(-1.0*Tmp1)*(180.0/PI);
        
		Tmp1 = vDOT(N23,N13)/(L1*L2);
		ang3 = acos(Tmp1)*(180.0/PI);
        
        fprintf(ahandle,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",jj,t,L1,L2,L3,ang1,ang2,ang3,R1/AU,R2/AU,R3/AU);
        fprintf(s1handle,"%d\t%g\t%g\t%g\t%g\t%g\t%g\n",jj,t,th1*(180.0/PI),R1/AU,s1[0],s1[1],s1[2]);
        fprintf(s2handle,"%d\t%g\t%g\t%g\t%g\t%g\t%g\n",jj,t,th2*(180.0/PI),R2/AU,s2[0],s2[1],s2[2]);
        fprintf(s3handle,"%d\t%g\t%g\t%g\t%g\t%g\t%g\n",jj,t,th3*(180.0/PI),R3/AU,s3[0],s3[1],s3[2]);
        
        
        
        //-----------------------------------------------------------
        //---------- GET ORBITAL PARAMETERS -------------------------
        //-----------------------------------------------------------
        
        //--- Look for ASCENDING NODES ------
        if ((s1old[2] < 0.0) && (s1[2] > 0.0))
           GetNode(1,t,delT,Rearth,ecc,Torb,kg,beta1,&Tanode1,&THanode1,&Ranode1,&Omega1);

        if ((s2old[2] < 0.0) && (s2[2] > 0.0))
           GetNode(1,t,delT,Rearth,ecc,Torb,kg,beta2,&Tanode2,&THanode2,&Ranode2,&Omega2);

        if ((s3old[2] < 0.0) && (s3[2] > 0.0))
           GetNode(1,t,delT,Rearth,ecc,Torb,kg,beta3,&Tanode3,&THanode3,&Ranode3,&Omega3);
        
        
        //--- Look for DESCENDING NODES ------
        if ((s1old[2] > 0.0) && (s1[2] < 0.0))
           GetNode(0,t,delT,Rearth,ecc,Torb,kg,beta1,&Tdnode1,&THdnode1,&Rdnode1,&Tmp1);

        if ((s2old[2] > 0.0) && (s2[2] < 0.0))
           GetNode(0,t,delT,Rearth,ecc,Torb,kg,beta2,&Tdnode2,&THdnode2,&Rdnode2,&Tmp1);

        if ((s3old[2] > 0.0) && (s3[2] < 0.0))
           GetNode(0,t,delT,Rearth,ecc,Torb,kg,beta3,&Tdnode3,&THdnode3,&Rdnode3,&Tmp1);
        
        
        //--- Look for APOCENTER ------
        if ((drdt1old > 0.0) && (drdt1 < 0.0))
           GetExtremum(0,t,delT,Rearth,ecc,Torb,kg,beta1,R1,&Tapo1,&THapo1,&Rapo1);

        if ((drdt2old > 0.0) && (drdt2 < 0.0))
           GetExtremum(0,t,delT,Rearth,ecc,Torb,kg,beta2,R2,&Tapo2,&THapo2,&Rapo2);

        if ((drdt3old > 0.0) && (drdt3 < 0.0))
           GetExtremum(0,t,delT,Rearth,ecc,Torb,kg,beta3,R3,&Tapo3,&THapo3,&Rapo3);
        
        
        //--- Look for PERICENTER ------
        if ((drdt1old < 0.0) && (drdt1 > 0.0))
           GetExtremum(1,t,delT,Rearth,ecc,Torb,kg,beta1,R1,&Tperi1,&THperi1,&Rperi1);

        if ((drdt2old < 0.0) && (drdt2 > 0.0))
           GetExtremum(1,t,delT,Rearth,ecc,Torb,kg,beta2,R2,&Tperi2,&THperi2,&Rperi2);

        if ((drdt3old < 0.0) && (drdt3 > 0.0))
           GetExtremum(1,t,delT,Rearth,ecc,Torb,kg,beta3,R3,&Tperi3,&THperi3,&Rperi3);

        
     }  //END LOOP: jj = nPTS
     
     fclose(ahandle);
     fclose(s1handle);
     fclose(s2handle);
     fclose(s3handle);
     
     
     //------------------------------------------------------------
     // Compute ARG OF PERIHELION from ASC NODE and PERICENTER     
     //------------------------------------------------------------
     
     w1 = THperi1 - THanode1;          //angle between ASC NODE and PERICENTER
     if (w1 < 0.0) w1 += (2.0*PI);   //make angle positive valued           
     
     w2 = THperi2 - THanode2;          //angle between ASC NODE and PERICENTER
     if (w2 < 0.0) w2 += (2.0*PI);   //make angle positive valued           
     
     w3 = THperi3 - THanode3;          //angle between ASC NODE and PERICENTER
     if (w3 < 0.0) w3 += (2.0*PI);   //make angle positive valued           
     
     
     
     
     //------------------------------------------------------------------------
     //--- Output Orbital Parameters ------------------------------------------
     //------------------------------------------------------------------------
     
     
     printf("SPACECRAFT 1 --------------------------------------------------\n");
     printf("apocenter = %g AU   TH = %g deg   t = %g s\n",Rapo1/AU,THapo1*(180.0/PI),Tapo1);
     printf("pericenter = %g AU   TH = %g deg   t = %g s\n",Rperi1/AU,THperi1*(180.0/PI),Tperi1);
     printf("ASC NODE   R = %g AU   TH = %g deg   t = %g s\n",Ranode1/AU,THanode1*(180.0/PI),Tanode1);
     printf("DSC NODE   R = %g AU   TH = %g deg   t = %g s\n",Rdnode1/AU,THdnode1*(180.0/PI),Tdnode1);
     printf("w = %g     OMEGA = %g deg\n\n",w1*180.0/PI,Omega1*180.0/PI);
     
     printf("SPACECRAFT 2 --------------------------------------------------\n");
     printf("apocenter = %g AU   TH = %g deg   t = %g s\n",Rapo2/AU,THapo2*(180.0/PI),Tapo2);
     printf("pericenter = %g AU   TH = %g deg   t = %g s\n",Rperi2/AU,THperi2*(180.0/PI),Tperi2);
     printf("ASC NODE   R = %g AU   TH = %g deg   t = %g s\n",Ranode2/AU,THanode2*(180.0/PI),Tanode2);
     printf("DSC NODE   R = %g AU   TH = %g deg   t = %g s\n",Rdnode2/AU,THdnode2*(180.0/PI),Tdnode2);
     printf("w = %g     OMEGA = %g deg\n\n",w2*180.0/PI,Omega2*180.0/PI);
     
     printf("SPACECRAFT 3 --------------------------------------------------\n");
     printf("apocenter = %g AU   TH = %g deg   t = %g s\n",Rapo3/AU,THapo3*(180.0/PI),Tapo3);
     printf("pericenter = %g AU   TH = %g deg   t = %g s\n",Rperi3/AU,THperi3*(180.0/PI),Tperi3);
     printf("ASC NODE   R = %g AU   TH = %g deg   t = %g s\n",Ranode3/AU,THanode3*(180.0/PI),Tanode3);
     printf("DSC NODE   R = %g AU   TH = %g deg   t = %g s\n",Rdnode3/AU,THdnode3*(180.0/PI),Tdnode3);
     printf("w = %g     OMEGA = %g deg\n\n",w3*180.0/PI,Omega3*180.0/PI);
     
     
   
   
   printf("\n\nALL DONE!  Whoo hoo!\n");
   
   return 0;
}



//==================================================================================
//----------------------------------------------------------------------------------
//==================================================================================


// ***************************************************************************
//                                                                         
//FUNCTION:  GetNode                                                       
//                                                                         
//Bisection routine for finding time and anomaly of an orbital node.       
//                                                                         
// ***************************************************************************

void GetNode(int kASCEND, double time, double delT, double Rearth, double ecc,
             double Torb, double kg, double beta, double *nodeTime, double *nodeAnomaly,
             double *nodeR, double *Omega)
{
  //============= VARIABLE DECLARATIONS =============
  
  double tLo, tHi, t1, th;
  double sc1[3], balpha;
  double eps, u;
  
  int Nb, bFLAG;
  
  double Tmp1;
  
  //=============== START ROUTINE HERE ===============
  
  tLo = time - delT;   //set the bracket values we are searching between
  tHi = time;          //based on the input values from the orbit       
  
  Nb = 0;              //initialize bisection counter
  eps = 1.0e-6;        //bisection tolerance value   
  
  bFLAG = 0;           //set flag to make sure we bisect at least once
  
  while (bFLAG == 0)
  {
     t1 = tLo + (tHi - tLo)/2.0;
     
     //compute spacecraft data at time = t1
      balpha = 2.0*PI*t1/Torb + kg;   //RCP alpha
      GetSpacecraftLocation(sc1,balpha,beta,ecc,Rearth);
      Tmp1 = (balpha - beta);
      th = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     u = sc1[2];            //spacecraft z value at t1
     

     //--- Descending node case ------
     if (kASCEND == 0)
     {
        if (u > 0.0) tLo = t1; //if under t, move lower bracket up 
        
        if (u < 0.0) tHi = t1; //if over t, move upper bracket down
     }
     
     
     //--- Ascending node case ------
     if (kASCEND == 1)
     {
        if (u < 0.0) tLo = t1; //if under t, move lower bracket up 
        
        if (u > 0.0) tHi = t1; //if over t, move upper bracket down
     }
     
     
     if ( fabs(u) <= eps ) bFLAG = 1;     //if reached tolerance, kick out of bisect
     
     if (Nb > 10000)                      //Terminate for a bisection error         
      {
         printf("ERROR in GetNode(): Bisections exceed 10000\n");
         break;
      }
      
      Nb++;
  
  } //<--- end WHILE <---
  
  //return values
  *Omega = atan2(sc1[1],sc1[0]);
  *nodeR = sqrt(vDOT(sc1,sc1));
  *nodeTime = t1;
  *nodeAnomaly = th;
  
  return;
}


// ***************************************************************************
//                                                                         
//FUNCTION:  GetExtremum                                                   
//                                                                         
//Bisection routine for finding time and anomaly of an orbital node.       
//                                                                         
// ***************************************************************************

void GetExtremum(int kPERI, double time, double delT, double Rearth, double ecc,
             double Torb, double kg, double beta, double Rvin, double *extrmTime,
             double *extrmAnomaly, double *extrmR)
{
  //============= VARIABLE DECLARATIONS =============
  
  double tLo, tHi, t1, th;
  double sc1[3], balpha;
  double R1, Rvold, drdtv;
  double eps, u;
  
  int Nb, bFLAG;
  
  double Tmp1;
  
  //=============== START ROUTINE HERE ===============
  
  tLo = time - delT;   //set the bracket values we are searching between
  tHi = time;          //based on the input values from the orbit       
  
  R1 = Rvin;           //Current s/c barycentric radius
  
  Nb = 0;              //initialize bisection counter
  eps = 1.0e-6;        //bisection tolerance value   
  
  bFLAG = 0;           //set flag to make sure we bisect at least once
  
  while (bFLAG == 0)
  {
     Rvold = R1;
     
     t1 = tLo + (tHi - tLo)/2.0;
     
     //compute spacecraft data at time = t1
      balpha = 2.0*PI*t1/Torb + kg;   //RCP alpha
      GetSpacecraftLocation(sc1,balpha,beta,ecc,Rearth);
      Tmp1 = (balpha - beta);
      th = Tmp1 + 2.0*ecc*sin(Tmp1) + (5.0/2.0)*ecc*ecc*cos(Tmp1)*sin(Tmp1);
     
     R1 = sqrt(vDOT(sc1,sc1));
     drdtv = (R1 - Rvold)/delT;

     u = drdtv;            //set u = to drdt; when drdt = 0 then at turning point
     

     //--- Apocenter case ------
     if (kPERI == 0)
     {
        if (u < 0.0) tLo = t1; //if under t, move lower bracket up 
        
        if (u > 0.0) tHi = t1; //if over t, move upper bracket down
     }
     
     
     //--- Pericenter case ------
     if (kPERI == 1)
     {
        if (u > 0.0) tLo = t1; //if under t, move lower bracket up 
        
        if (u < 0.0) tHi = t1; //if over t, move upper bracket down
     }
     
     
     if ( fabs(u) <= eps ) bFLAG = 1;     //if reached tolerance, kick out of bisect
     
     if (Nb > 10000)                      //Terminate for a bisection error         
      {
         printf("ERROR in GetNode(): Bisections exceed 10000\n");
         break;
      }
      
      Nb++;
  
  } //<--- end WHILE <---
  
  //return values
  *extrmR = sqrt(vDOT(sc1,sc1));
  *extrmTime = t1;
  *extrmAnomaly = th;
  
  return;
}


// ***************************************************************************
//                                                                         
//FUNCTION:  GetSpacecraftLocation                                         
//                                                                         
//Rubbo-Cornish-Poujade model for s/c location in LISA constellation.      
//                                                                         
// ***************************************************************************

void GetSpacecraftLocation(double baryctr[],double a, double b, double e, double R)
{
  //============= VARIABLE DECLARATIONS =============
  
  double Tmp1, Tmp2;
  
  //=============== START ROUTINE HERE ===============
  
  Tmp1 = R*cos(a) + (e/2.0)*R*(cos(2.0*a - b) - 3.0*cos(b));
  Tmp2 = (e*e/8.0)*R*(3.0*cos(3.0*a - 2.0*b) - 10.0*cos(a) - 5.0*cos(a - 2.0*b));
  baryctr[0] = Tmp1 + Tmp2;
  
  Tmp1 = R*sin(a) + (e/2.0)*R*(sin(2.0*a - b) - 3.0*sin(b));
  Tmp2 = (e*e/8.0)*R*(3.0*sin(3.0*a - 2.0*b) - 10.0*sin(a) + 5.0*sin(a - 2.0*b));
  baryctr[1] = Tmp1 + Tmp2;
  
  Tmp1 = -1.0*sqrt(3.0)*e*R*cos(a - b);
  Tmp2 = sqrt(3.0)*e*e*R*(cos(a - b)*cos(a - b) + 2.0*sin(a - b)*sin(a - b));
  baryctr[2] = Tmp1 + Tmp2;
  
}




// ***************************************************************************
//                                                                         
//FUNCTION:  vDOT                                                          
//TYPE    :  double                                                        
//                                                                         
//This routine computes the 3D dot product of 2 vectors.                   
//                                                                         
//PARAMETERS                                                               
//-------------------------------------------------------------------------
//                                                                         
// ***************************************************************************

double vDOT(double v1[], double v2[])
{
  //============= VARIABLE DECLARATIONS =============
  
  int i;
  
  double ans;
  
  //=============== START ROUTINE HERE ===============
  
  ans = 0.0;
  
  for(i = 0; i <= 2; i++)
     ans += v1[i]*v2[i];
  
  return ans;
}


// ***************************************************************************
//                                                                         
//FUNCTION:  Die                                                           
//TYPE    :  int                                                           
//RETURNS :  0 (void)                                                      
//MODIFIED:  15 February 2002                                              
//                                                                         
//'Die' error messges then gracefully exits the program.                   
//                                                                         
// ***************************************************************************

void Die(char routine[], char errorMsg[])
{
 //=============== VARIABLES ===============

 //=============== START ROUTINE HERE ===============
   
   printf("\n\n\n");
   printf("DYING in routine %s\n",routine);
   printf("PROBLEM: %s\n\n",errorMsg);
   
   fflush(stdin);
   printf("Hit a key -- Hasta la vista, baby.");
   getchar();
   
   exit;
   
   return;

} //FUNCTION: DIE DIE DIE
