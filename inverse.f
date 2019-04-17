!****************************************** Inverse **************************************************
!
!  This code is used to calculate the inverse matrix of square matirx COEFF.
!  Version 1.0  2017. 11 .12  created by Dan Xue.
!
!*****************************************************************************************************
!     Input  
!       COEFF : The matrix needs to inverse 
!		  N	: The dimension of matirx COEFF
!     Output
!        ACOF : The inverse matrix of matirx COEFF 
!*****************************************************************************************************

      SUBROUTINE inverse(COEFF,ACOF,N)
      INTEGER N
      REAL COEFF(N,N),ACOF(N,N)   
      INTEGER IP(N)      
      REAL P          
      INTEGER I0,R      
	ACOF = COEFF
      EPS  = 0.0000000001

      DO K = 1,N
          P  = 0.
          I0 = K
          IP(K) = K
          
          DO I = K,N
              IF(ABS(ACOF(I,K)).GT.ABS(P))THEN
                  P  = ACOF(I,K)
                  I0 = I
                  IP(K) = I
              ENDIF
          ENDDO
          
          IF(ABS(P).LE.EPS)THEN
              WRITE(*,*) 'DET=0'
              stop 
          ENDIF
          
          IF(I0.NE.K)THEN
              DO J  = 1,N
                  S = ACOF(K,J)
                  ACOF(K,J)  = ACOF(I0,J)
                  ACOF(I0,J) = S
              ENDDO
          ENDIF
         
          ACOF(K,K) = 1./P
          
          DO I = 1,N
              IF(I.NE.K)THEN
                  ACOF(I,K) = -ACOF(I,K)*ACOF(K,K)
                  DO J = 1,N
                      IF(J.NE.K)THEN
                          ACOF(I,J) = ACOF(I,J)+ACOF(I,K)*ACOF(K,J)
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
         
          DO J = 1,N
              IF(J.NE.K)THEN
                  ACOF(K,J) = ACOF(K,K)*ACOF(K,J)
              ENDIF
          ENDDO
       ENDDO

       DO K = N-1,1,-1
          R = IP(K)
          IF(R.NE.K)THEN
              DO I  = 1,N
                  S = ACOF(I,R)
                  ACOF(I,R) = ACOF(I,K)
                  ACOF(I,K) = S
              ENDDO
          ENDIF
      ENDDO

      END