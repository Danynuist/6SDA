!************************* Azimuth-dependent single scattering approximation ************************
!
!  The code is used to calculate the  bidirectional reflectance of single scattering approximation
!  which is azimuth-dependent. In addition, the delta-adjustment is not used.
!   Version 1.0  2017. 02 .13  created by Feng Zhang
!
!****************************************************************************************************
!     Input:
!         Opt_Dep(K) : Optical  depth   of the k-th sub-layer 
!             SSA(K) : Single scattering albedo  of the k-th sub-layer    
!       Cof_Leg(K,M) : The expansion coefficient of phase function of the k-th sub-layer
!                       Cof_Leg(K,1)=SSA1=3*g, Cof_Leg(K,2)=SSA2,....
!             U0_Sol : Cosine of the  solar zenith angle
!              U_Sat : Cosine of the light exiting zenith angle
!            Del_Phi : Del_Phi=varph_{0}-varph,  the soalr azimuth - the light exiting azimuth
!                Lay : Number of sub-layers 
!     Output:
!            Phi_SSA : Bidirectional reflectance at U_Sat and Del_Phi      Liou book 3.4.11bIN put:
!                      Intensity I(0, U_Sat, Del_Phi) =Phi_SSA * U0_Sol * F0 / PI
!  Variable :
!        Opt_Dep_Acc : Accumulate delta-adjustment optical depth
!            Pha_Fun : Azimuth-dependent phase function 
!****************************************************************************************************
	SUBROUTINE  single_sca_azi(Opt_Dep, SSA, Cof_Leg, Num_Leg, 
     &U0_Sol, U_Sat, Del_Phi, Lay, Phi_SSA)
      IMPLICIT REAL (A-H,O-Z)
	PARAMETER  (PI = 3.1415927)
      REAL Cof_Leg(Lay,Num_Leg), Opt_Dep(Lay), SSA(Lay), OTAU(Lay),
     &     P(Num_Leg+1)
      INTEGER j 

      D9      = 1.0 / 13.0
      Co_ta   = -1 * U_Sat * U0_Sol + (1-U_Sat ** 2 ) ** 0.5     
     &        * (1-U0_Sol ** 2) ** 0.5 * COS (Del_Phi)
         
	Opt_Dep_Acc = 0.0
      Phi_SSA     = 0.0

      DO 100  K   = 1, Lay
	   F        = D9 * Cof_Leg(K,6)
         OTAU(K)  = Opt_Dep(K) * (1.0 - SSA(K) * F)
         
	   P(1)     = 1
	   P(2)     = Co_ta

	   DO j     = 3,257
            P(j)  = ((2*j-3)*P(j-1)*Co_ta-(j-2)*P(j-2))/(j-1)
 	   END DO

         Pha_Fun  = 1

	   DO i     =1,256
	      Pha_Fun = Pha_Fun + p(i+1) * Cof_Leg(K,i)
	   END DO

	 
         Phi_SSA  = Phi_SSA 
     &            + SSA(K) / (4 * U_Sat + 4 * U0_Sol) * Pha_Fun
     &            * (1 - EXP(- OTAU (k) * (1.0/U_Sat + 1.0/U0_Sol)))
     &            * EXP(- Opt_Dep_Acc/U_Sat ) * EXP(-Opt_Dep_Acc/U0_Sol)
     &            / (1.0 - SSA(K) * F)   
         Opt_Dep_Acc = OTAU (k) + Opt_Dep_Acc
        
     
  100 CONTINUE

      END SUBROUTINE 

!*********************** Azimuthally averaged single scattering approximation **********************
! 
!  The code is used to calculate the  azimuthally averaged bidirectional reflectance of single scattering approximation
!  In addition, the delta-adjustment is used. Positive exponent Exp (k*Opt_Dep) is eliminated in the code
!   
!  Version 1.0  2017. 02 .13  created by Feng Zhang
!****************************************************************************************************
!     Input:
!         Opt_Dep(K) : Optical  depth   of the k-th sub-layer 
!             SSA(K) : Single scattering albedo  of the k-th sub-layer    
!       Cof_Leg(K,M) : Expansion coefficient of phase function of the k-th sub-layer
!                       Cof_Leg(K,1)=SSA1=3*g, Cof_Leg(K,2)=SSA2,....
!             U0_Sol : Cosine of the  solar zenith angle
!              U_Sat : Cosine of the light exiting zenith angle
!                Lay : The number of sub-layers 
!     Output:
!              U_SSA : Azimuthally averaged bidirectional reflectance at U_Sat     
!                      Intensity I(0, U_Sat) =U_SSA * U0_Sol * F0 / PI
!     Variable:       
!        Opt_Dep_Acc : Accumulate delta-adjustment optical depth
!            Pha_Fun : Azimuth-independent phase function 
!****************************************************************************************************
      SUBROUTINE single_sca_azi_ave(Opt_Dep, SSA, Cof_Leg, U0_Sol, 
     &                              U_Sat, Lay, U_SSA)
      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (M = 256)
	PARAMETER  (PI = 3.1415927)
      REAL Cof_Leg(Lay,M), Opt_Dep(Lay), SSA(Lay), OTAU(Lay)
        
	D9          = 1.0 / 13.0
	Opt_Dep_Acc = 0.0
	U_SSA       = 0.0
      DO 100  K     = 1, Lay
	   F          = D9 * Cof_Leg(K,6) 
         OM         = SSA(K) * (1.0 - F) / (1.0 - SSA(K) * F)
         OTAU(K)    = Opt_Dep(K) * (1.0 - SSA(K) * F)

         OMEGA1     = OM * (Cof_Leg(K,1) - 3.0 * F) / (1.0 - F)
         OMEGA2     = OM * (Cof_Leg(K,2) - 5.0 * F) / (1.0 - F)
         OMEGA3     = OM * (Cof_Leg(K,3) - 7.0 * F) / (1.0 - F)
	   OMEGA4     = OM * (Cof_Leg(K,4) - 9.0 * F) / (1.0 - F)
	   OMEGA5     = OM * (Cof_Leg(K,5) - 11.0 * F) / (1.0 - F)
         Pha_Fun    = 1 - OMEGA1  * U0_Sol * U_Sat + OMEGA2 * 0.5    
     &              * (3 * U0_Sol**2 - 1.0) * 0.5 * (3 * U_Sat**2 - 1.0)
     &              - OMEGA3 * 0.5 * ( 5 * U0_Sol  ** 3 -  3 *  U0_Sol)
     &              * 0.5 * ( 5 * U_Sat  ** 3 -  3 *  U_Sat)
     &              + OMEGA4 * 0.125 * ( 35 * U0_Sol**4 -30 * U0_Sol**2
     &              + 3 ) * 0.125 * ( 35 * U_Sat**4 -30 * U_Sat**2 + 3) 
     &              - OMEGA5 * 0.125 * ( 63 * U0_Sol **5 -70 * U0_Sol**3
     &              + 15 * U0_Sol ) * 0.125 * ( 63 * U_Sat**5 
     &              - 70 * U_Sat**3 + 15 * U_Sat ) 
        
         U_SSA      = U_SSA 
     &              + OM  / (4 * U_Sat + 4 * U0_Sol) * Pha_Fun
     &              * (1 - EXP(- OTAU (k) * (1/U_Sat + 1/U0_Sol)))
     &              * EXP(- Opt_Dep_Acc/ U_Sat) 
     &              * EXP(- Opt_Dep_Acc/ U0_Sol )    
         Opt_Dep_Acc = OTAU(k) + Opt_Dep_Acc

    
  100 CONTINUE

      END SUBROUTINE 
