!************************************  Adding method ************************************************ 
!
!  The code is used to calculate the reflectance and I0 to I5 which used in azimuth-dependent 
!  intensity calculation at the top of the atmosphere based on reflection and tranmission of each
!  sub-layer. And transmission calculation is included.
!
!****************************************************************************************************
!     Input
!               Lay : The number of sub-layer 
!           Flx_Sol : The solar flux
!            U0_Sol : Cosine of solar zenith angle  
!              Rbar : Matrices for diffuse reflection  
!	         Tbar : Matrices for diffuse transmission
!	          Ru0 : Matrix for direct reflection
!	          Tu0 : Matrix for direct transmission
!	       Dir_Tr : Direct transmission Dir_Tr = EXP(-Opt_Dep/U0_Sol)
!     Output
!         TOA_6S_I0 : I0 at the top of the atmosphere in 6SDA
!         TOA_6S_I1 : I1 at the top of the atmosphere in 6SDA
!         TOA_6S_I2 : I2 at the top of the atmosphere in 6SDA
!         TOA_6S_I3 : I3 at the top of the atmosphere in 6SDA
!         TOA_6S_I4 : I4 at the top of the atmosphere in 6SDA
!         TOA_6S_I5 : I5 at the top of the atmosphere in 6SDA
!            Re_TOA : Reflection at the top of the atmosphere ( you can output reflection in this code 
!                     as you wish, but here we don't output )
!            Tr_SFC : Transmission at the bottom of the atmosphere ( you can output transmission in this
!                     code as you wish, but here we don't output )
!*****************************************************************************************************
	SUBROUTINE add6s(Lay, Lev, Flx_Sol, U0_Sol,
     &                 Rbar, Tbar, Ru0, Tu0, Dir_Tr,
     &                 TOA_6S_I0, TOA_6S_I1, TOA_6S_I2, 
     &                 TOA_6S_I3, TOA_6S_I4, TOA_6S_I5,Re_TOA,Tr_SFC) 
      IMPLICIT REAL (A-H,O-Z),
     &INTEGER (I-N)
      PARAMETER (Lay1=1, Lev1=1, PI=3.1415927)								   
	REAL Ru0(Lev,3,1), Tu0(Lev,3,1), Dir_Tr(Lev),
     &     Rbar(Lev,3,3), Tbar(Lev,3,3)
	REAL Rbar1(Lev,3,3), Tbar1(Lev,3,3), 
     &     Rbars1(Lev,3,3), Tbars1(Lev,3,3), 
     &     RbarN(Lev,3,3), TbarN(Lev,3,3), 
     &     Ru01(Lev,3,1), Tu01(Lev,3,1),
     &     Ru0N(Lev,3,1), Tu0N(Lev,3,1), Dir_Tr1(Lev)
	REAL Flxu(Lev), Flxd(Lev)
      REAL EYE(3,3), TMM(3,3), TMR(3,3), 
     &     BMR(3,3), BRT(3,3),
     &     Uu0(3,1), Du0(3,1), COEFF(3,3)
     	DATA EYE/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/  

	F      = U0_Sol * Flx_Sol
	Albedo = 0.0
      Trans  = 0.0

	Rbar(Lev,1,1) = Albedo
      Rbar(Lev,2,1) = - 0.25 * Albedo
	Rbar(Lev,3,1) = 0.125 * Albedo

      Rbar(Lev,1:3,2:3) = 0.0
      Tbar(Lev,1:3,1:3) = Trans

      Ru0(Lev,1,1) = Albedo
	Ru0(Lev,2,1) = - 0.25 * Albedo
	Ru0(Lev,3,1) = 0.125 * Albedo

	Tu0(Lev,1:3,1) = Trans

	Dir_Tr(Lev) = Trans
	                  
	Rbar1(Lay1,1:3,1:3)  = Rbar(Lay1,1:3,1:3)
	Tbar1(Lay1,1:3,1:3)  = Tbar(Lay1,1:3,1:3) 
	Rbars1(Lay1,1:3,1:3) = Rbar(Lay1,1:3,1:3)
	Tbars1(Lay1,1:3,1:3) = Tbar(Lay1,1:3,1:3)
	 
      Ru01(Lay1,:,:) = Ru0(Lay1,:,:)
	Tu01(Lay1,:,:) = Tu0(Lay1,:,:)
	Dir_Tr1(Lay1)    = Dir_Tr(Lay1)

      DO 100 k = 2,Lev
	   COEFF = EYE - MATMUL(Rbar(k,:,:),Rbars1(k-1,:,:))
         CALL inverse(COEFF,TMM,3)
         BMR = MATMUL(TMM,Tbar(k,:,:))

         Dir_Tr1(k) = Dir_Tr1(k-1) * Dir_Tr(k) 
         Uu0(:,1) = MATMUL(TMM,Ru0(k,:,1) * Dir_Tr1(k-1)              
     &            + MATMUL(Rbar(k,:,:),Tu01(k-1,:,1)))
	   Du0(:,1) = Tu01(k-1,:,1) + MATMUL(Rbars1(k-1,:,:),Uu0(:,1)) 	  
         Tu01(k,:,1) = Tu0(k,:,1) * Dir_Tr1(k-1)
     &               + MATMUL(Tbar(k,:,:),Du0(:,1))

         BRT(:,:) = MATMUL(Tbar(k,:,:),Rbars1(k-1,:,:))
	   Rbars1(k,:,:) = Rbar(k,:,:) + MATMUL(BRT(:,:),BMR)     
 100	CONTINUE

         Ru0N(Lev,1:3,1)    = Ru0(Lev,1:3,1)
     	   Tu0N(Lev,1:3,1)    = Tu0(Lev,1:3,1)
	   RbarN(Lev,1:3,1:3) = Rbar(Lev,1:3,1:3)
	   TbarN(Lev,1:3,1:3) = Tbar(Lev,1:3,1:3) 

      DO 200 k = Lay,Lay1,-1
	   COEFF = EYE - MATMUL(RbarN(k+1,:,:),Rbar(k,:,:))
         CALL inverse(COEFF,TMM,3)
            
         Uu0(:,1) = MATMUL(TMM, Ru0N(k+1,:,1) * Dir_Tr(k)         
     &            + MATMUL(RbarN(k+1,:,:),Tu0(k,:,1)))
	   Du0(:,1) = Tu0(k,:,1) + MATMUL(Rbar(k,:,:),Uu0(:,1)) 
	   Ru0N(k,:,1) = Ru0(k,:,1) + MATMUL(Tbar(k,:,:),Uu0(:,1))        

	   COEFF = EYE - MATMUL(Rbar(k,:,:),RbarN(k+1,:,:))
         CALL inverse(COEFF,TMM,3)      
         TMR = MATMUL(Tbar(k,:,:),MATMUL(RbarN(k+1,:,:),TMM(:,:)))
         RbarN(k,:,:) = Rbar(k,:,:) + MATMUL(TMR,Tbar(k,:,:))
200   CONTINUE

         Flxu(Lev1) = Ru0N(1,1,1) * F
	   Flxd(Lev1) = F

	   F1 = Flxu(1)
	   F2 = Ru0N(1,2,1) * F
	   F3 = Ru0N(1,3,1) * F

      DO 300 K = Lev1,Lev-1
	   COEFF = EYE - MATMUL(RbarN(k+1,:,:),Rbars1(k,:,:))
         CALL inverse(COEFF,TMM,3)

	   Uu0(:,1) = MATMUL(TMM, Ru0N(k+1,:,1) * Dir_Tr1(k)          
     &            + MATMUL(RbarN(k+1,:,:),Tu01(k,:,1)))
	   Du0(:,1) = Tu01(k,:,1) + MATMUL(Rbars1(k,:,:),Uu0(:,1)) 
         Flxu(k+1) = Uu0(1,1) * F
         Flxd(k+1) = Du0(1,1) * F + Dir_Tr1(k) * F
300   CONTINUE
	    
      Re_TOA = Flxu(1) / F     
      Tr_SFC = Flxd(Lev) / F   

	TOA_6S_I4 = (32./63)*(F1/(4*PI)- (F1-F2)/(5*PI)-(F1-8*F3)/(14*PI)) 
      TOA_6S_I0 = (8./5)  *((F1-F2)  / (4*PI)+105./128*TOA_6S_I4)              
	TOA_6S_I2 = (16./35)*((F1-8*F3)/ (4*PI)+84. /16 *TOA_6S_I4)                
	TOA_6S_I1 = F1 / (4*PI)                                               
	TOA_6S_I3 = F2 / (4*PI)                                               
	TOA_6S_I5 = F3 / (4*PI)
                                               
	END SUBROUTINE  

