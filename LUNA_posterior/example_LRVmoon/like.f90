MODULE like

use params
use transitmod
use radialmod

implicit none
      
contains
      
!=======================================================================

SUBROUTINE slikelihood(R,slhood)
         
	implicit none
      
	REAL(8), DIMENSION(nest_nPar) :: R
        REAL(8) :: slhood, loglike
	INTEGER :: i

        ! Scaling of the parameters from hypercube
	DO i = 1, sdim
          IF( Rflag(i) .EQ. 0 .OR. Rflag(i) .EQ. 2 ) THEN ! Modes 0 and 2
            ! Uniform: Rmax = max, Rmin = min
	    R(i) = Rmin(i) + (Rmax(i)-Rmin(i))*R(i)
          ELSE IF( Rflag(i) .EQ. 3 ) THEN ! Mode 3
            ! Gaussian: Rmax = mean, Rmin = stdev
            R(i) = Rmax(i) + roottwo*Rmin(i)*inverf(-1.0D0+2.0D0*R(i))
          ELSE IF( Rflag(i) .EQ. 4 ) THEN ! Mode 4
            ! Jeffrey's: Rmax = max, Rmin = min
            R(i) = ( Rmax(i)**R(i) )*( Rmin(i)**(1.0D0-R(i)) )
          ELSE IF( Rflag(i) .EQ. 5 ) THEN ! Mode 5
            ! Modified Jefrey's: Rmax = max, Rmin = inflection point
            R(i) = -( Rmin(i)**(1.0D0-R(i)) )*( Rmin(i)**R(i) - ( Rmin(i)+Rmax(i) )**R(i) )
          END IF
	END DO

        ! Call transit to get chi^2
        call models(R,loglike,0)

	slhood = loglike
        !!write(*,*) 'R(1) = ',R(1),' yielding chi2 = ',chi2

END SUBROUTINE slikelihood
      
!=======================================================================

!! ======================================================================
SUBROUTINE models(Rvec,loglike,showpri)

!implicit none

 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 INTEGER :: showpri, showrv
 REAL(8) :: loglike, loglikeR
 INTEGER :: i, j
 REAL(8) :: jitter
 REAL(8) :: chi2RV
 REAL(8), DIMENSION(nplen) :: resP
 REAL(8) :: loglikeP, loglikeP_lc, loglikeP_sc
 INTEGER :: nplen_lc, nplen_sc
 REAL(8), DIMENSION(nslen) :: resS
 REAL(8) :: loglikeS, loglikeS_lc, loglikeS_sc
 INTEGER :: nslen_lc, nslen_sc
 REAL(8), DIMENSION(nqlen) :: resQ
 REAL(8) :: loglikeQ, loglikeQ_lc, loglikeQ_sc
 INTEGER :: nqlen_lc, nqlen_sc
 REAL(8), DIMENSION(ntlen) :: resT
 REAL(8) :: loglikeT, loglikeT_lc, loglikeT_sc
 INTEGER :: ntlen_lc, ntlen_sc
 REAL(8), DIMENSION(nulen) :: resU
 REAL(8) :: loglikeU, loglikeU_lc, loglikeU_sc
 INTEGER :: nulen_lc, nulen_sc

 ! === Call transit for Kepler data ===

 call transit(Rvec,nplen,resP,showpri,&
              tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,&
              NresamP,integ_bigP,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamP .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeP = 0.0D0
     DO i=1,nplen
       loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
     END DO
     loglikeP = -0.5D0*nplen*LogTwoPi - 0.5D0*loglikeP
   ELSE
     ! == LC ==
     loglikeP_lc = 0.0D0
     nplen_lc = 0
     DO i=1,nplen
       IF( sigfp(i) .LE. flick ) THEN !Take only the LC data
         loglikeP_lc = loglikeP_lc + resP(i)**2*sigfpwei(i) + logsigfp(i)
         nplen_lc = nplen_lc + 1
       END IF
     END DO
     loglikeP_lc = -0.5D0*nplen_lc*LogTwoPi - 0.5D0*loglikeP_lc
     ! == SC ==
     loglikeP_sc = 0.0D0
     nplen_sc = 0
     DO i=1,nplen
       IF( sigfp(i) .GT. flick ) THEN !Take only the SC data
         loglikeP_sc = loglikeP_sc + resP(i)**2*sigfpwei(i) + logsigfp(i)
         nplen_sc = nplen_sc + 1
       END IF
     END DO
     loglikeP_sc = -0.5D0*nplen_sc*LogTwoPi - 0.5D0*loglikeP_sc
     ! == BOTH ==
     loglikeP = loglikeP_lc + loglikeP_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeP = 1.0D0 
 END IF

 
 ! === Call transit for TESS data ===
 
 call transit(Rvec,nslen,resS,showpri,&
              ts,fs,sigfs,epochs_seq,fswei,sigfswei,&
              NresamS,integ_bigS,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamS .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeS = 0.0D0
     DO i=1,nslen
       loglikeS = loglikeS + resS(i)**2*sigfswei(i) + logsigfs(i)
     END DO
     loglikeS = -0.5D0*nslen*LogTwoPi - 0.5D0*loglikeS
   ELSE
     ! == LC ==
     loglikeS_lc = 0.0D0
     nslen_lc = 0
     DO i=1,nslen
       IF( sigfs(i) .LE. flick ) THEN !Take only the LC data
         loglikeS_lc = loglikeS_lc + resS(i)**2*sigfswei(i) + logsigfs(i)
         nslen_lc = nslen_lc + 1
       END IF
     END DO
     loglikeS_lc = -0.5D0*nslen_lc*LogTwoPi - 0.5D0*loglikeS_lc
     ! == SC ==
     loglikeS_sc = 0.0D0
     nslen_sc = 0
     DO i=1,nslen
       IF( sigfs(i) .GT. flick ) THEN !Take only the SC data
         loglikeS_sc = loglikeS_sc + resS(i)**2*sigfswei(i) + logsigfs(i)
         nslen_sc = nslen_sc + 1
       END IF
     END DO
     loglikeS_sc = -0.5D0*nslen_sc*LogTwoPi - 0.5D0*loglikeS_sc
     ! == BOTH ==
     loglikeS = loglikeS_lc + loglikeS_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeS = 1.0D0 
 END IF
 
 ! === Call transit for BARON data ===

 call transit(Rvec,nqlen,resQ,showpri,&
              tq,fq,sigfq,epochq_seq,fqwei,sigfqwei,&
              NresamQ,integ_bigQ,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamQ .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeQ = 0.0D0
     DO i=1,nqlen
       loglikeQ = loglikeQ + resQ(i)**2*sigfqwei(i) + logsigfq(i)
     END DO
     loglikeQ = -0.5D0*nqlen*LogTwoPi - 0.5D0*loglikeQ
   ELSE
     ! == LC ==
     loglikeQ_lc = 0.0D0
     nqlen_lc = 0
     DO i=1,nqlen
       IF( sigfq(i) .LE. flick ) THEN !Take only the LC data
         loglikeQ_lc = loglikeQ_lc + resQ(i)**2*sigfqwei(i) + logsigfq(i)
         nqlen_lc = nqlen_lc + 1
       END IF
     END DO
     loglikeQ_lc = -0.5D0*nqlen_lc*LogTwoPi - 0.5D0*loglikeQ_lc
     ! == SC ==
     loglikeQ_sc = 0.0D0
     nqlen_sc = 0
     DO i=1,nqlen
       IF( sigfq(i) .GT. flick ) THEN !Take only the SC data
         loglikeQ_sc = loglikeQ_sc + resQ(i)**2*sigfqwei(i) + logsigfq(i)
         nqlen_sc = nqlen_sc + 1
       END IF
     END DO
     loglikeQ_sc = -0.5D0*nqlen_sc*LogTwoPi - 0.5D0*loglikeQ_sc
     ! == BOTH ==
     loglikeQ = loglikeQ_lc + loglikeQ_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeQ = 1.0D0 
 END IF
 
 ! === Call transit to LCO transit ===

 call transit(Rvec,ntlen,resT,showpri,&
              tt,ft,sigft,epocht_seq,ftwei,sigftwei,&
              NresamT,integ_bigT,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamT .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeT = 0.0D0
     DO i=1,ntlen
       loglikeT = loglikeT + resT(i)**2*sigftwei(i) + logsigft(i)
     END DO
     loglikeT = -0.5D0*ntlen*LogTwoPi - 0.5D0*loglikeT
   ELSE
     ! == LC ==
     loglikeT_lc = 0.0D0
     ntlen_lc = 0
     DO i=1,ntlen
       IF( sigft(i) .LE. flick ) THEN !Take only the LC data
         loglikeT_lc = loglikeT_lc + resT(i)**2*sigftwei(i) + logsigft(i)
         ntlen_lc = ntlen_lc + 1
       END IF
     END DO
     loglikeT_lc = -0.5D0*ntlen_lc*LogTwoPi - 0.5D0*loglikeT_lc
     ! == SC ==
     loglikeT_sc = 0.0D0
     ntlen_sc = 0
     DO i=1,ntlen
       IF( sigft(i) .GT. flick ) THEN !Take only the SC data
         loglikeT_sc = loglikeT_sc + resT(i)**2*sigftwei(i) + logsigft(i)
         ntlen_sc = ntlen_sc + 1
       END IF
     END DO
     loglikeT_sc = -0.5D0*ntlen_sc*LogTwoPi - 0.5D0*loglikeT_sc
     ! == BOTH ==
     loglikeT = loglikeT_lc + loglikeT_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeT = 1.0D0 
 END IF
 
 ! === Call transit to LCO transit ===

 call transit(Rvec,nulen,resU,showpri,&
              tu,fu,sigfu,epochu_seq,fuwei,sigfuwei,&
              NresamU,integ_bigU,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamU .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeU = 0.0D0
     DO i=1,nulen
       loglikeU = loglikeU + resU(i)**2*sigfuwei(i) + logsigfu(i)
     END DO
     loglikeU = -0.5D0*nulen*LogTwoPi - 0.5D0*loglikeU
   ELSE
     ! == LC ==
     loglikeU_lc = 0.0D0
     nulen_lc = 0
     DO i=1,nulen
       IF( sigfu(i) .LE. flick ) THEN !Take only the LC data
         loglikeU_lc = loglikeU_lc + resU(i)**2*sigfuwei(i) + logsigfu(i)
         nulen_lc = nulen_lc + 1
       END IF
     END DO
     loglikeU_lc = -0.5D0*nulen_lc*LogTwoPi - 0.5D0*loglikeU_lc
     ! == SC ==
     loglikeU_sc = 0.0D0
     nulen_sc = 0
     DO i=1,nulen
       IF( sigfu(i) .GT. flick ) THEN !Take only the SC data
         loglikeU_sc = loglikeU_sc + resU(i)**2*sigfuwei(i) + logsigfu(i)
         nulen_sc = nulen_sc + 1
       END IF
     END DO
     loglikeU_sc = -0.5D0*nulen_sc*LogTwoPi - 0.5D0*loglikeU_sc
     ! == BOTH ==
     loglikeU = loglikeU_lc + loglikeU_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeU = 1.0D0 
 END IF

 ! === Call radial ===

 IF ( rvflag .EQ. 1 ) THEN

 jitter = 0.0D0
 call radial(Rvec,chi2RV)

 ! Compute log likelihood under simple conditions of global noise terms
 loglikeR = chi2RV
 DO i=1,nrlen
   loglikeR = loglikeR + DLOG( sigrv(i)**2 + jitter**2 )
 END DO
 loglikeR = -0.5D0*nrlen*LogTwoPi - 0.5D0*loglikeR
 
 ELSE
   loglikeR = 0.0D0
 END IF

 ! === Sum up all loglikes (product of likes) ===
 loglike = loglikeP + loglikeR + loglikeS + loglikeQ + loglikeT + loglikeU
 !write(*,*) loglikeP,loglikeS,loglikeQ,loglikeT,loglikeU

END SUBROUTINE models
! ======================================================================

! ======================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ======================================================================

END MODULE like

