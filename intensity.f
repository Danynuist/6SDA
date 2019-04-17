!************************************* Intensity calculation *****************************************
!
!  The code is used to calculate the  bidirectional reflectance of 6SDA+SingleScattering approximation
!  which is azimuth-dependent based on result of adding method.  
!
!*****************************************************************************************************
!     Input
!         TOA_6S_I0 : I0 at the top of the atmosphere in 6SDA
!         TOA_6S_I1 : I1 at the top of the atmosphere in 6SDA
!         TOA_6S_I2 : I2 at the top of the atmosphere in 6SDA
!         TOA_6S_I3 : I3 at the top of the atmosphere in 6SDA
!         TOA_6S_I4 : I4 at the top of the atmosphere in 6SDA
!         TOA_6S_I5 : I5 at the top of the atmosphere in 6SDA
!               Lay : The number of sub-layer 
!           Num_Leg : Number of Legendre polynomials for expanding the phase function
!           Flx_Sol : The solar flux
!            U0_Sol : Cosine of solar zenith angle
!             U_Sat : Cosine of satellite's zenith angle
!           Phi_Sol : Solar azmuith angle (choose Pi means angle equal to 180)
!           Phi_Sat : Satellite azmuith angle (choose Pi means angle equal to 180)
!        Opt_Dep(K) : Optical  depth   of the k-th sub-layer 
!            SSA(K) : Single scattering albedo  of the k-th sub-layer
!      Cof_Leg(K,M) : Cofficients of the M-th Legendre polynomials
!                     For example, in HG phase function Cof_Leg(K,1)=SSA1=3*Asy_Fac(K),
!                     Cof_Leg(K,2)=SSA2=5*Asy_Fac(K)*Asy_Fac(K),....
!    Output
!            U_6S_I : bidirectional reflectance of 6SDA+SingleScattering approximation
!                     which is azimuth-dependent 
!*****************************************************************************************************
      SUBROUTINE intensity(TOA_6S_I0, TOA_6S_I1, TOA_6S_I2, 
     &                     TOA_6S_I3, TOA_6S_I4, TOA_6S_I5,
     &                     Lay, Num_Leg, Flx_Sol, U0_Sol, U_Sat, 
     &                     Phi_Sol, Phi_Sat, Opt_Dep, SSA, Cof_Leg,
     &                     U_6S_I)
      IMPLICIT REAL (A-H,O-Z),
     &INTEGER (I-N)
      PARAMETER (PI=3.1415927)					   
	REAL Opt_Dep(Lay), SSA(Lay), Cof_Leg(Lay,1:Num_Leg), U_6S_I
	

	TAU0  = SUM(Opt_Dep(:))              ! total optical depth
 	TAU0S = SUM(Opt_Dep(:) * SSA(:))     ! total scattering optical depth
    SSAS  = TAU0S/TAU0                   ! total Single scattering albedo
	
	PP0u = 1
	PP1u = U_Sat
	PP2u = 0.5*(3*U_Sat**2-1)
	PP3u = 0.5*(5*U_Sat**3-3*U_Sat)
	PP4u = 0.125*(35*U_Sat**4-30*U_Sat**2+3)
	PP5u = 0.125*(63*U_Sat**5-70*U_Sat**3+15*U_Sat)
	PP6u = 0.0625*(231*U_Sat**6-315*U_Sat**4+105*U_Sat**2-5)
	PP7u = (13*U_Sat*PP6u-6*PP5u)/7.

    A2 = 0.1996 + 0.0775 * EXP( -0.5842 * TAU0S )
	B2 = -0.3258 - 3.6672 * EXP( -1.0781 * TAU0S )
	C2 = 1.5586 + 2.0655 * EXP( -1.0068 * TAU0S )
	D2 = 0.3137 + 0.3078 * EXP( -1.6669 * TAU0S )

	AI = TOA_6S_I0+3*PP1u*TOA_6S_I1+5*PP2u*TOA_6S_I2+9*PP4u*TOA_6S_I4
     &   + (7*PP3u-SQRT(7.0)*1.6*A2*(PP4u-B2*PP5u)-SQRT(7.0)*A2
     &   * ((1-U_Sat)**4))*TOA_6S_I3
     &   + (11*PP5u-SQRT(11.0)*(16./7)*C2*(PP6u-D2*PP7u)
     &   + SQRT(11.0)*C2*((1-U_Sat)**6))*TOA_6S_I5
	 	  
      
	Del_Phi = Phi_Sol-Phi_Sat

	Phi_SSA_I = 0
	U_SSA_I   = 0
	Phi_6S_I  = 0
	 
      CALL SINGLE_SCA_AZI(Opt_Dep, SSA, Cof_Leg, Num_Leg, 
     &                    U0_Sol, U_Sat, Del_Phi, Lay, Phi_SSA)

      CALL SINGLE_SCA_AZI_AVE(Opt_Dep, SSA, Cof_Leg, U0_Sol, 
     &                        U_Sat, Lay, U_SSA)


	Phi_SSA_I = Phi_SSA * U0_Sol * Flx_Sol / PI
	U_SSA_I   = U_SSA   * U0_Sol * Flx_Sol / PI
      

    AI1 = Phi_SSA_I - U_SSA_I + AI
		
		
      IF ( TAU0 >= 1. ) THEN 	
        A1=(1.13210807803445-1.44785466368487*SSAS)
     &     *(1-exp(-1.18546534815896*TAU0))*exp(-5.10258975067384*U_Sat)
     &     *sin(-2.42168344004181*U_Sat)
        B1=(-17.3311616102198+22.8904124676226*SSAS)
     &     *(1-exp(-1.26466476568627*TAU0))*exp(-5.51147313593757*U_Sat)
     &     *sin(0.0956664267093702*U_Sat)
        C1=(-0.198249994180936+0.273685349949774*SSAS)
     &     *(1-exp(-1.72152600300149*TAU0))*exp(-6.20917598029264*U_Sat)
     &     *sin(1.90165492156087*U_Sat)
        D1=(-0.0421096158792793+0.0596389591331587*SSAS)
     &     *(1-exp(-2.17787697666655*TAU0))*exp(-6.80988538046632*U_Sat)
     &     *sin(2.95470932954989*U_Sat)
        E1=(0.198553092873448-0.17239421200655*SSAS)
     &     *(1-exp(-18891.7215032446*TAU0))*exp(-8.71205917473957*U_Sat)
     &     *sin(9.82447000568429*U_Sat)  
	 
	 ELSE IF (TAU0<1.) THEN
         A1=(-1.6807474459889+2.34116968451252*SSAS)
     &      *exp(-7.6589634199901*U_Sat)*sin(-0.115239508363325*U_Sat)
     &      *(9.10534206844836*TAU0+32.571061567654*TAU0**2
     &      -51.2437639453175*TAU0**3+22.0522948702216*TAU0**4
     &      -0.122246901909939)
	 
         B1=(-7.50491479890676+10.6745197625413*SSAS)
     &      *exp(-7.13701770515571*U_Sat)*sin(0.0678798912939125*U_Sat)
     &      *(-2.16486527332583*TAU0-7.26238747868187*TAU0**2
     &      +11.709020335856*TAU0**3-5.10100509313766*TAU0**4
     &      +0.0287055115576268)

        C1=(-1.4860828376762+2.16282258370169*SSAS)
     &      *exp(-8.07958328001287*U_Sat)*sin(-0.0239731488331223*U_Sat)
     &      *(12.6894332189626*TAU0+24.9961416468932*TAU0**2
     &      -49.9480718633997*TAU0**3+23.6369227475064*TAU0**4
     &      -0.162141146016741)
	 
         D1=(-0.321453371896571+0.474028411853402*SSAS)
     &      *exp(-8.6228882228047*U_Sat)*sin(1.96064490030737*U_Sat)
     &      *(-0.347599764629956*TAU0-0.332293662053404*TAU0**2
     &      +0.935455746483804*TAU0**3-0.483909530157814*TAU0**4
     &      +0.00429334517406901)
	 
          E1=(0.00407312763990981+0.00823747806213116*SSAS)
     &      *exp(-12.701051630485*U_Sat)*sin(-10.6977082544266*U_Sat)
     &      *(32.335245895412*TAU0-53.7816475671741*TAU0**2
     &      +31.8672035426173*TAU0**3-5.06723151786867*TAU0**4
     &      -0.122387955720022)
	 ENDIF
	 
      U_6S_I  = AI1+A1*cos(Phi_Sat)+B1*cos(2*Phi_Sat)+
     &           C1*cos(3*Phi_Sat)+D1*cos(4*Phi_Sat)+E1
     

	END SUBROUTINE
