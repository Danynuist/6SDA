!*********************** 6SDA combined with single scattering approximation ************************
!
!  The code is used to calculate the  bidirectional reflectance of 6SDA+SingleScattering approximation
!  which is azimuth-dependent.
!  Version 1.0  2017. 11 .12  created by Feng Zhang, Dan Xue and Yi-Ning Shi
!
!****************************************************************************************************
!     Input parameter (detail in the input parameter section)
!              Lay : The number of sub-layer 
!          Num_Leg : Number of Legendre polynomials for expanding the phase function
!          Flx_Sol : The solar flux
!           U0_Sol : Cosine of solar zenith angle
!            U_Sat : Cosine of satellite's zenith angle
!          Phi_Sol : Solar azmuith angle (choose Pi means angle equal to 180)
!          Phi_Sat : Satellite azmuith angle (choose Pi means angle equal to 180)
!      Get_Pha_Fun : Choose for phase function, 1 is phase function; 2 is Cloud C phase function
!
!     Input profile (detail in the input profile section)
!       Opt_Dep(K) : Optical depth of the k-th sub-layer 
!           SSA(K) : Single scattering albedo  of the k-th sub-layer
!       Asy_Fac(K) : Scattering asymmetry factor of the k-th sub-layer used in the HG phase function      
!     Cof_Leg(K,M) : Cofficients of the M-th Legendre polynomials
!                    For example, in HG phase function Cof_Leg(K,1)=SSA1=3*Asy_Fac(K), 
!                    Cof_Leg(K,2)=SSA2=5*Asy_Fac(K)*Asy_Fac(K),....
! 
!     Output 
!           U_6S_I : bidirectional reflectance of 6SDA+SingleScattering approximation which is azimuth-dependent
!****************************************************************************************************   
      PROGRAM six spheric stream adding
      IMPLICIT REAL (A-H,O-Z),
     &INTEGER (I-N)
!************************************* Input parameter **********************************************
      PARAMETER (Lay=1, Lay1=1, Lev1=1, Lev=Lay+1, Num_Leg=256,
     &           PI=3.1415927, Flx_Sol=1.0, U0_Sol=0.5, U_Sat=0.5,
     &           Phi_Sol=0, Phi_Sat=0, Get_Pha_Fun=2)
!      'Num_Leg' is Number of Legendre polynomials for expanding the phase function
!      'Flx_Sol' is solar flux
!      'U0_Sol' is cosine of solar zenith angle
!      'U_Sat' is cosine of satellite's zenith angle
!      'Phi_Sol' is solar azmuith angle (choose PI means angle equal to 180)
!      'Phi_Sat' is satellite azmuith angle (choose PI means angle equal to 180)
!      'Get_Pha_Fun' is choose for phase function, 1 is phase function; 2 is Cloud C phase function
!****************************************************************************************************
	REAL Opt_Dep(Lay), SSA(Lay), Asy_Fac(Lay), Cof_Leg(Lay,1:Num_Leg), 
     &     WWW(Num_Leg)
      REAL Ru0(Lev,3,1), Tu0(Lev,3,1), Dir_Tr(Lev),
     &     Rbar(Lev,3,3), Tbar(Lev,3,3),tau(32),u(21),phi(5),tau1(53)

!*************************************** Input profile ***********************************************
     	DO i = Lay1, Lay
         Opt_Dep(i) =10  ! optical depth for each sub-layer
    	   SSA(i) = 0.99! single-scattering albedo for each sub-layer
		   
    	   IF ( Get_Pha_Fun == 1 ) THEN
!----------------- HG approximation -------------
             Asy_Fac(i) = 0.8   ! scattering asymmetry factor for each sub-layer 
    	       DO j = 1, Num_Leg
	           Cof_Leg(i, j) = (2. * j + 1) * Asy_Fac(i) ** (FLOAT(j))
     	       END DO	 
!------------------------------------------------
         ELSE IF ( Get_Pha_Fun == 2 ) THEN
!----------------- Cloud C ----------------------
    	       DO j = 1, Num_Leg
    	           OPEN(19, file = "Cof_Leg_cl.txt")
    	           READ(19, *) WWW(j)
    	           Cof_Leg(i,j) = (2. * j + 1) * WWW(j)
			   ! cofficients of each Legendre polynomials   
     	       END DO
	       CLOSE(19)
!------------------------------------------------
          END IF
      END DO
!****************************************************************************************************

      DO i = Lay1, Lay
	CALL sixspherical(Opt_Dep(i), SSA(i), Cof_Leg(i,1:6), U0_Sol,
     &                 Rbar(i,:,:), Tbar(i,:,:), Ru0(i,:,1), Tu0(i,:,1),
     &                 Dir_Tr(i))
	END DO    
	!6SDA for each sub-layer

      CALL add6s(Lay, Lev, Flx_Sol, U0_Sol,
     &           Rbar, Tbar, Ru0, Tu0, Dir_Tr,
     &           TOA_6S_I0, TOA_6S_I1, TOA_6S_I2, 
     &           TOA_6S_I3, TOA_6S_I4, TOA_6S_I5, Re_TOA,Tr_SFC)
      !Adding method based on reflection and transmission of all sub-layer

	CALL intensity(TOA_6S_I0, TOA_6S_I1, TOA_6S_I2, 
     &               TOA_6S_I3, TOA_6S_I4, TOA_6S_I5,
     &               Lay, Num_Leg, Flx_Sol, U0_Sol, U_Sat, 
     &               Phi_Sol, Phi_Sat, Opt_Dep, SSA, Cof_Leg,
     &               U_6S_I)
	 
	!calculating the bidirectional reflectance at the top of atmosphere, which is azimuth-dependent 

!****************************************************************************************************
      END PROGRAM
