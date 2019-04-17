!*********************************** 6SDA for each sub-layer ****************************************
!
!  The code is used to calculate the reflection and transmission for each layer by six-streams spherical 
!  harmonic expansion approximation.
!
!****************************************************************************************************
!     Input:  
!        Opt_Dep : Optical  depth 
!            SSA : Single scattering albedo
!        Cof_Leg : Cofficients of the Legendre polynomials
!	    U0_Sol : Cosine of solar zenith angle
!
!     Output:
!           Rbar : Matrices for diffuse reflection  
!	      Tbar : Matrices for diffuse transmission
!	       Ru0 : Matrix for direct reflection
!	       Tu0 : Matrix for direct transmission
!	    Dir_Tr : Direct transmission Dir_Tr = EXP(-Opt_Dep/U0_Sol)
!****************************************************************************************************

      SUBROUTINE  sixspherical( Opt_Dep, SSA, Cof_Leg, U0_Sol,
     &                          Rbar, Tbar, Ru0, Tu0, Dir_Tr )
	IMPLICIT REAL (A-H,O-Z),
     &INTEGER (I-N)
      PARAMETER (M=6, PI=3.1415927)
      REAL Cof_Leg(M),OMEGA(M)
	REAL Rbar(1:3,1:3), Tbar(1:3,1:3), Ru0(1:3,1), Tu0(1:3,1)
	 
	REAL A(0:M-1), B0(0:M-1), COEFF(1:M,1:M), INV(6,6), ALAM(1:3),
     &     P(1:3), Q(1:3), R(1:3), Y(1:3), Z(1:3), G(1:6),
     &     SU(1:3), SD(1:3), TU(1:3), TD(1:3), WU(1:3), WD(1:3), TR(1:3)

!************************************** Delta function adjustment ***********************************
      F   = Cof_Leg(6) / 13.
      DO I=1,6
         OMEGA(I) = Cof_Leg(I)
      END DO

      OM  = SSA * (1.-F) / (1.-SSA*F)  !SSA after adjustment      
      TAU = Opt_Dep * (1.-SSA*F)       !Opt_Dep after adjustment

	A(0)=1.-OM  
      DO L = 1,M-1
         OMEGA(L)=(OMEGA(L)-(2.*L+1.)*F)/(1.-F) !Cof_Leg after adjustment
         A(L)=(2.*L+1.)-OM*OMEGA(L)       
      END DO  
!****************************************************************************************************

      B0(0) =  OM/ 4.                         
      B0(1) = -OM* OMEGA(1)* U0_Sol / 4.          
      B0(2) =  OM* OMEGA(2)* (3. *U0_Sol**2-1.)/8.
      B0(3) = -OM* OMEGA(3)* (5. *U0_Sol**3-3. *U0_Sol) /8.
      B0(4) =  OM* OMEGA(4)* (35.*U0_Sol**4-30.*U0_Sol**2+3) /32.
      B0(5) = -OM* OMEGA(5)* (63.*U0_Sol**5-70.*U0_Sol**3+15*U0_Sol)/32.
   
      UU = A(0)*A(1) + 4./9.*A(0)*A(3) + 64./225.*A(0)*A(5)
     &   + 1./9.*A(2)*A(3) + 16./225.*A(2)*A(5) + 1./25.*A(4)*A(5)
      VV = 1./9.*A(0)*A(1)*A(2)*A(3)  + 16./225.*A(0)*A(1)*A(2)*A(5)
     &   + 1./25.*A(0)*A(1)*A(4)*A(5) + 4./225.*A(0)*A(3)*A(4)*A(5)
     &   + 1./225.*A(2)*A(3)*A(4)*A(5)
      WW = 1./225.*A(0)*A(1)*A(2)*A(3)*A(4)*A(5)
      
      AP = UU**2-3.*VV    
      BP = -UU*VV+9.*WW							
      THETA = ACOS((-2.*AP*UU-3.*BP)/(2.*SQRT(AP**3)))

      ALAM(1) = SQRT((UU - 2. * SQRT(AP) * COS(THETA/3.)) / 3.)			
      ALAM(2) = SQRT((UU + SQRT(AP) * (COS(THETA/3.)
     &        + SQRT(3.) * SIN(THETA/3.))) / 3.)	   
      ALAM(3) = SQRT((UU + SQRT(AP) * (COS(THETA/3.)
     &        - SQRT(3.) * SIN(THETA/3.))) / 3.)		
      
	AA = 1./U0_Sol 
      
      Dir_Tr   = EXP(-TAU*AA)
      TR(1) = EXP(-ALAM(1)*TAU)
      TR(2) = EXP(-ALAM(2)*TAU)
      TR(3) = EXP(-ALAM(3)*TAU)

      DELTA = -225. * (AA**6 - UU*AA**4 + VV*AA**2 - WW)
      DELTA = 1./DELTA
      
      ETA0=((-225.*B0(1)+150.*B0(3)-120.*B0(5))*AA**5-(-225.*A(1)*B0(0)
     &    -100.*A(3)*B0(0)+50.*A(3)*B0(2)-64.*A(5)*B0(0)+32.*A(5)*B0(2)
     &    -24.*A(5)*B0(4))*AA**4-(-25.*A(2)*A(3)*B0(1)-16.*A(2)*A(5)
     &    *B0(1)-9.*A(4)*A(5)*B0(1)+6.*A(4)*A(5)*B0(3))*AA**3-(25.*A(1)
     &    *A(2)*A(3)*B0(0)+16.*A(1)*A(2)*A(5)*B0(0)+9.*A(1)*A(4)*A(5)
     &    *B0(0)+4.*A(3)*A(4)*A(5)*B0(0)-2.*A(3)*A(4)*A(5)*B0(2))*AA**2
     &    -A(2)*A(3)*A(4)*A(5)*B0(1)*AA+A(1)*A(2)*A(3)*A(4)*A(5)*B0(0))
     &    *DELTA
      
      ETA1=(-225.*B0(0)*AA**5+(225.*A(0)*B0(1)-150.*A(0)*B0(3)+120.*A(0)
     &    *B0(5))*AA**4+(50.*A(0)*A(3)*B0(2)+25.*A(2)*A(3)*B0(0)+32.
     &    *A(0)*A(5)*B0(2)+16.*A(2)*A(5)*B0(0)-24.*A(0)*A(5)*B0(4)+9.
     &    *A(4)*A(5)*B0(0))*AA**3+(-25.*A(0)*A(2)*A(3)*B0(1)-16.*A(0)
     &    *A(2)*A(5)*B0(1)-9.*A(0)*A(4)*A(5)*B0(1)+6.*A(0)*A(4)*A(5)
     &    *B0(3))*AA**2+(-2.*A(0)*A(3)*A(4)*A(5)*B0(2)-A(2)*A(3)*A(4)
     &    *A(5)*B0(0))*AA+A(0)*A(2)*A(3)*A(4)*A(5)*B0(1))*DELTA
     
      ETA2=(-(75.*B0(3)-60.*B0(5))*AA**5-(50.*A(3)*B0(0)-25.*A(3)*B0(2)
     &    +32.*A(5)*B0(0)-16.*A(5)*B0(2)+12.*A(5)*B0(4))*AA**4-(-75.
     &    *A(0)*A(1)*B0(3)-50.*A(0)*A(3)*B0(1)+60.*A(0)*A(1)*B0(5)-32.
     &    *A(0)*A(5)*B0(1)-3.*A(4)*A(5)*B0(3))*AA**3-(25.*A(0)*A(1)*A(3)
     &    *B0(2)+16.*A(0)*A(1)*A(5)*B0(2)-12.*A(0)*A(1)*A(5)*B0(4)-2.
     &    *A(3)*A(4)*A(5)*B0(0)+A(3)*A(4)*A(5)*B0(2))*AA**2-(3.*A(0)
     &    *A(1)*A(4)*A(5)*B0(3)+2.*A(0)*A(3)*A(4)*A(5)*B0(1))*AA+A(0)
     &    *A(1)*A(3)*A(4)*A(5)*B0(2))*DELTA
     
      ETA3=((150.*B0(0)-75.*B0(2))*AA**5+(-150.*A(0)*B0(1)+100.*A(0)
     &    *B0(3)-80.*A(0)*B0(5)+25.*A(2)*B0(3)-20.*A(2)*B0(5))*AA**4
     &    +(75.*A(0)*A(1)*B0(2)+16.*A(0)*A(5)*B0(4)-6.*A(4)*A(5)*B0(0)
     &    +4.*A(2)*A(5)*B0(4)+3.*A(4)*A(5)*B0(2))*AA**3+(-25.*A(0)*A(1)
     &    *A(2)*B0(3)+20.*A(0)*A(1)*A(2)*B0(5)+6.*A(0)*A(4)*A(5)*B0(1)
     &    -4.*A(0)*A(4)*A(5)*B0(3)-A(2)*A(4)*A(5)*B0(3))*AA**2+(-4.*A(0)
     &    *A(1)*A(2)*A(5)*B0(4)-3.*A(0)*A(1)*A(4)*A(5)*B0(2))*AA+A(0)
     &    *A(1)*A(2)*A(4)*A(5)*B0(3))*DELTA
      
      ETA4=(-45.*B0(5)*AA**5+(24.*A(5)*B0(0)-12.*A(5)*B0(2)+9.*A(5)
     &    *B0(4))*AA**4+(45.*A(0)*A(1)*B0(5)-24.*A(0)*A(5)*B0(1)+20.
     &    *A(0)*A(3)*B0(5)+16.*A(0)*A(5)*B0(3)+5.*A(2)*A(3)*B0(5)+4.
     &    *A(2)*A(5)*B0(3))*AA**3+(12.*A(0)*A(1)*A(5)*B0(2)-9.*A(0)*A(1)
     &    *A(5)*B0(4)-4.*A(0)*A(3)*A(5)*B0(4)-A(2)*A(3)*A(5)*B0(4))
     &    *AA**2+(-5.*A(0)*A(1)*A(2)*A(3)*B0(5)-4.*A(0)*A(1)*A(2)*A(5)
     &    *B0(3))*AA+A(0)*A(1)*A(2)*A(3)*A(5)*B0(4))*DELTA
     
      ETA5=(-(120.*B0(0)-60.*B0(2)+45.*B0(4))*AA**5-(-120.*A(0)*B0(1)
     &    +80.*A(0)*B0(3)-64.*A(0)*B0(5)+20.*A(2)*B0(3)-16.*A(2)*B0(5)
     &    -9.*A(4)*B0(5))*AA**4-(60.*A(0)*A(1)*B0(2)-45.*A(0)*A(1)*B0(4)
     &    -20.*A(0)*A(3)*B0(4)-5.*A(2)*A(3)*B0(4))*AA**3-(-20.*A(0)*A(1)
     &    *A(2)*B0(3)+16.*A(0)*A(1)*A(2)*B0(5)+9.*A(0)*A(1)*A(4)*B0(5)
     &    +4.*A(0)*A(3)*A(4)*B0(5)+A(2)*A(3)*A(4)*B0(5))*AA**2-5.*A(0)
     &    *A(1)*A(2)*A(3)*B0(4)*AA+A(0)*A(1)*A(2)*A(3)*A(4)*B0(5))*DELTA   
      
      DO i = 1,3
          P(i) = -1. / a(0)*(a(2)*a(3)*a(4)*a(5)/(120.*ALAM(i)**3)
     &         - 5.*a(2)*a(3)/(24.*ALAM(i))-2.*a(2)*a(5)/(15.*ALAM(i))
     &         - 3.*a(4)*a(5)/(40.*ALAM(i))+15.*ALAM(i)/8.)
	    Q(i) = a(2)*a(3)*a(4)*a(5)/(120.*ALAM(i)**4)
     &         - 5.*a(2)*a(3)/(24*ALAM(i)**2)-2.*a(2)*a(5)
     &         / (15.*ALAM(i)**2)-3.*a(4)*a(5)/(40*ALAM(i)**2)+15./8.
	    R(i) = -a(3)*a(4)*a(5)/(60.*ALAM(i)**3)+5.*a(3)/(12.*ALAM(i))
     &         + 4.*a(5)/(15.*ALAM(i))
	    Y(i) = a(4)*a(5)/(20.*ALAM(i)**2)-5./4.
	    Z(i) = -a(5)/(5.*ALAM(i))
      END DO
      
      H1 = -( 1./2.*ETA0-ETA1+5./8.*ETA2-3./16.*ETA4)
      H2 = -(-1./8.*ETA0+5./8.*ETA2-ETA3+81./128.*ETA4)
      H3 = -(1./16.*ETA0-25./128.*ETA2+81./128.*ETA4-1.*ETA5)
      H4 = -( 1./2.*ETA0+ETA1+5./8.*ETA2-3./16.*ETA4)*Dir_Tr
      H5 = -(-1./8.*ETA0+5./8.*ETA2 +ETA3+81./128.*ETA4)*Dir_Tr
      H6 = -(1./16.*ETA0-25./128.*ETA2 +81./128.*ETA4+1.*ETA5)*Dir_Tr
                   
      DO i=1,3
          SU(i)= 1./2.*P(i)+Q(i)+5./8.*R(i)-3./16.*Z(i)
          SD(i)= 1./2.*P(i)-Q(i)+5./8.*R(i)-3./16.*Z(i)
          TU(i)=-1./8.*P(i)+5./8.*R(i)+Y(i)+81./128.*Z(i)
          TD(i)=-1./8.*P(i)+5./8.*R(i)-Y(i)+81./128.*Z(i)
          WU(i)=1./16.*P(i)-25./128.*R(i)+81./128.*Z(i)+1
          WD(i)=1./16.*P(i)-25./128.*R(i)+81./128.*Z(i)-1
      END DO
      
      DO i  = 1,6
          N = MOD(i,2)
          IF (N.NE.0) THEN
          L = INT(i/2) + MOD(i,2)
          COEFF(1,i) = SD(L)
          COEFF(2,i) = TD(L)
          COEFF(3,i) = WD(L)
          COEFF(4,i) = SU(L) * TR(L)
          COEFF(5,i) = TU(L) * TR(L)
          COEFF(6,i) = WU(L) * TR(L)
          ELSE
          L=INT(i/2)
          COEFF(1,i) = SU(L) * TR(L)
          COEFF(2,i) = TU(L) * TR(L)
          COEFF(3,i) = WU(L) * TR(L)
          COEFF(4,i) = SD(L)
          COEFF(5,i) = TD(L)
          COEFF(6,i) = WD(L)
          END IF
      END DO

      CALL inverse(COEFF,INV,M)
      
      DO i=1,6
          G(i) = INV(i,1)*H1 + INV(i,2)*H2 + INV(i,3)*H3
     &         + INV(i,4)*H4 + INV(i,5)*H5 + INV(i,6)*H6					 
      END DO 
  
      TEU0 = P(1)*(G(1)+G(2)*TR(1)) + P(2)*(G(3)+G(4)*TR(2))
     &     + P(3)*(G(5)+G(6)*TR(3)) + ETA0
      TEU1 = Q(1)*(G(1)-G(2)*TR(1)) + Q(2)*(G(3)-G(4)*TR(2))
     &     + Q(3)*(G(5)-G(6)*TR(3)) + ETA1
      TEU2 = R(1)*(G(1)+G(2)*TR(1)) + R(2)*(G(3)+G(4)*TR(2))
     &     + R(3)*(G(5)+G(6)*TR(3)) + ETA2
      TEU3 = Y(1)*(G(1)-G(2)*TR(1)) + Y(2)*(G(3)-G(4)*TR(2))
     &     + Y(3)*(G(5)-G(6)*TR(3)) + ETA3
      TEU4 = Z(1)*(G(1)+G(2)*TR(1)) + Z(2)*(G(3)+G(4)*TR(2))
     &     + Z(3)*(G(5)+G(6)*TR(3)) + ETA4
	TEU5 = G(1)-G(2)*TR(1) + G(3)-G(4)*TR(2)
     &     + G(5)-G(6)*TR(3) + ETA5  
      

      Ru0(1,1) = 2.*( 1./2.*TEU0+TEU1+5./8.*TEU2-3./16.*TEU4)*AA
      Ru0(2,1) = 2.*(-1./8.*TEU0+5./8.*TEU2+TEU3+81./128.*TEU4)*AA
      Ru0(3,1) = 1./4.*( 1./2.*TEU0-25./16.*TEU2+81./16.*TEU4+8*TEU5)*AA
      
                   
      TED0 = P(1)*(G(1)*TR(1)+G(2)) + P(2)*(G(3)*TR(2)+G(4))
     &     + P(3)*(G(5)*TR(3)+G(6)) + ETA0*Dir_Tr 
      TED1 = Q(1)*(G(1)*TR(1)-G(2)) + Q(2)*(G(3)*TR(2)-G(4))
     &     + Q(3)*(G(5)*TR(3)-G(6)) + ETA1*Dir_Tr 
      TED2 = R(1)*(G(1)*TR(1)+G(2)) + R(2)*(G(3)*TR(2)+G(4))
     &     + R(3)*(G(5)*TR(3)+G(6)) + ETA2*Dir_Tr 
      TED3 = Y(1)*(G(1)*TR(1)-G(2)) + Y(2)*(G(3)*TR(2)-G(4))
     &     + Y(3)*(G(5)*TR(3)-G(6)) + ETA3*Dir_Tr 
      TED4 = Z(1)*(G(1)*TR(1)+G(2)) + Z(2)*(G(3)*TR(2)+G(4))
     &     + Z(3)*(G(5)*TR(3)+G(6)) + ETA4*Dir_Tr 
	TED5 = G(1)*TR(1)-G(2) + G(3)*TR(2)-G(4)
     &     + G(5)*TR(3)-G(6) + ETA5*Dir_Tr  
      

      Tu0(1,1)= 2.*( 1./2.*TED0-TED1+5./8.*TED2-3./16.*TED4)*AA
      Tu0(2,1)= 2.*(-1./8.*TED0+5./8.*TED2-TED3+81./128.*TED4)*AA
      Tu0(3,1)= 1./4.*( 1./2.*TED0-25./16.*TED2+81./16.*TED4-8*TED5)*AA
     
      H1 = 1
      H2 = 1
      H3 = 1
      
      DO j = 1,3
          DO i = 1,6
             G(i) = INV(i,j) * H1
          END DO           
          TEU0 = P(1)*(G(1)+G(2)*TR(1)) + P(2)*(G(3)+G(4)*TR(2))
     &         + P(3)*(G(5)+G(6)*TR(3))
          TEU1 = Q(1)*(G(1)-G(2)*TR(1)) + Q(2)*(G(3)-G(4)*TR(2))
     &         + Q(3)*(G(5)-G(6)*TR(3))
          TEU2 = R(1)*(G(1)+G(2)*TR(1)) + R(2)*(G(3)+G(4)*TR(2))
     &         + R(3)*(G(5)+G(6)*TR(3))
          TEU3 = Y(1)*(G(1)-G(2)*TR(1)) + Y(2)*(G(3)-G(4)*TR(2))
     &         + Y(3)*(G(5)-G(6)*TR(3))
          TEU4 = Z(1)*(G(1)+G(2)*TR(1)) + Z(2)*(G(3)+G(4)*TR(2))
     &         + Z(3)*(G(5)+G(6)*TR(3))
	    TEU5 = G(1)-G(2)*TR(1) + G(3)-G(4)*TR(2)
     &         + G(5)-G(6)*TR(3)  
      
          Rbar(1,j) = 1./2.*TEU0+TEU1+5./8.*TEU2-3./16.*TEU4
          Rbar(2,j) = -1./8.*TEU0+5./8.*TEU2+TEU3+81./128.*TEU4
          Rbar(3,j) = 1./16.*TEU0-25./128.*TEU2+81./128.*TEU4+1*TEU5	  
      
                    
          TED0 = P(1)*(G(1)*TR(1)+G(2)) + P(2)*(G(3)*TR(2)+G(4))
     &         + P(3)*(G(5)*TR(3)+G(6))
          TED1 = Q(1)*(G(1)*TR(1)-G(2)) + Q(2)*(G(3)*TR(2)-G(4))
     &         + Q(3)*(G(5)*TR(3)-G(6))
          TED2 = R(1)*(G(1)*TR(1)+G(2)) + R(2)*(G(3)*TR(2)+G(4))
     &         + R(3)*(G(5)*TR(3)+G(6)) 
          TED3 = Y(1)*(G(1)*TR(1)-G(2)) + Y(2)*(G(3)*TR(2)-G(4))
     &         + Y(3)*(G(5)*TR(3)-G(6))
          TED4 = Z(1)*(G(1)*TR(1)+G(2)) + Z(2)*(G(3)*TR(2)+G(4))
     &         + Z(3)*(G(5)*TR(3)+G(6))
	    TED5 = G(1)*TR(1)-G(2) + G(3)*TR(2)-G(4)
     &         + G(5)*TR(3)-G(6)

          Tbar(1,j) =  1./2.*TED0-TED1+5./8.*TED2-3./16.*TED4
          Tbar(2,j) = -1./8.*TED0+5./8.*TED2-TED3+81./128.*TED4
          Tbar(3,j) = 1./16.*TED0-25./128.*TED2+81./128.*TED4-1*TED5
      END DO

      END SUBROUTINE 
