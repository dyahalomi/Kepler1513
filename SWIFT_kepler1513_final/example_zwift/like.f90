MODULE like

use params
use transitmod

implicit none
      
CONTAINS
      
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

 ! Call models to get likelihood
 call models(R,loglike,0)

 slhood = loglike

END SUBROUTINE slikelihood
      
!=======================================================================

! ======================================================================
SUBROUTINE models(Rvec,loglike,showpri)

!implicit none

 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 INTEGER :: showpri, showrv
 ! REAL(8), DIMENSION(nplen) :: resP
 REAL(8), DIMENSION(MT) :: resP
 REAL(8) :: loglike, loglikeP, loglikeE
 INTEGER :: i, chaoticflag,badhillflag

 ! === Call transit to primary transit ===

!call transit(Rvec,nplen,resP,&
!              showpri,tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,&
!              NresamP,integ_bigP,&
!              .FALSE.)

 call transit(Rvec,resP)

 ! TTV model likelihood
 loglikeP = 0.0D0
 IF( sflag .EQ. 1 .OR. sflag .EQ. 2 ) THEN
   DO i=1,nt(1)
     loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
   END DO
   loglikeP = -0.5D0*nt(1)*LogTwoPi - 0.5D0*loglikeP
 END IF
 IF( sflag .EQ. 3 ) THEN
   DO i=1,2*nt(1)
     loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
   END DO
   loglikeP = -0.5D0*2.0d0*nt(1)*LogTwoPi - 0.5D0*loglikeP
 END IF
 
 ! Add on photoeccentric penalty
 loglikeE = photoeccentric(Rvec(5),Rvec(7),&
            -0.3167948578830856D0,0.1113429775302461D0,0.30632857778554495D0) ! mu,sigmaneg,sigmapos
			
 ! test if chaotic, inputs (P1,e1,w1,mu1,P2,e2,w2,mu2)
 chaoticflag = chaotic(Rvec(11),Rvec(5),Rvec(7),10.0D0**Rvec(1),&
                        Rvec(4),Rvec(6),Rvec(8),Rvec(2))
 ! test if hill unstable, inputs (P1,e1,mu1,P2,e2,mu2)
 badhillflag = hillunstable(Rvec(4),Rvec(6),Rvec(2),&
                            Rvec(11),Rvec(5),10.0D0**Rvec(1))
 
 ! === Sum up all loglikes (product of likes) ===
 IF( chaoticflag .EQ. 0 .AND. badhillflag .EQ. 0 ) THEN
   ! Add on eccentricity prior | transiting
   ! IMPORTANT - the constant of integration depends on alpha+beta
   ! but I am here just setting it as a constant for alpha=0.867,beta=3.03
   ! so you'll need to update if you change the prior shape
   loglike = loglikeE + loglikeP + &
             loglikeprior(Rvec(5),Rvec(7)) + &
			 loglikebetaprior(Rvec(6))
 ELSE
   loglike = -HUGE(1.0D0)
 END IF
 	
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

! ======================================================================
FUNCTION photoeccentric(eccb,varpi_deg,mu,sigmaneg,sigmapos)

 implicit none

 ! Inputs
 REAL(8) :: eccb, varpi_deg
 REAL(8) :: mu, sigmaneg, sigmapos
 ! Intermediate
 REAL(8) :: loglikePH, logPsib, omegab
 ! Constants
 REAL(8), PARAMETER :: twooverpi = 0.7978845608028654D0
 ! Output
 REAL(8) :: photoeccentric

 ! I *think* that omega_b = varpi_b - Omega_b and that Omega_b = 1.5pi by definition
 omegab = (varpi_deg - 270.0D0)*radian
 
 ! Now define log(Psi)
 logPsib = DLOG( ( (1.0D0+eccb*DSIN(omegab))**3 ) / ( (DSQRT(1.0D0-eccb**2))**3 ) )
 
 IF( logPsib .LE. mu ) THEN
   loglikePH = (logPsib - mu)**2 / (2.0D0*sigmaneg**2)
   loglikePH = twooverpi*DEXP( -loglikePH )
   loglikePH = DLOG( loglikePH/(sigmaneg+sigmapos) )
 ELSE
   loglikePH = (logPsib - mu)**2 / (2.0D0*sigmapos**2)
   loglikePH = twooverpi*DEXP( -loglikePH )
   loglikePH = DLOG( loglikePH/(sigmaneg+sigmapos) )
 END IF
 
 photoeccentric = loglikePH

END FUNCTION
! ======================================================================

! ======================================================================
FUNCTION loglikeprior(eccb,varpi_deg)

 implicit none

 ! Inputs
 REAL(8) :: eccb, varpi_deg
 ! Intermediate
 REAL(8) :: omegab
 ! Constants
 REAL(8), PARAMETER :: alpha0 = 0.867D0
 REAL(8), PARAMETER :: beta0 = 3.03D0
 REAL(8), PARAMETER :: intconst0 = 3.0365830138756844D0 ! needs to be updated by hand
 ! Output
 REAL(8) :: loglikeprior

 ! I *think* that omega_b = varpi_b - Omega_b and that Omega_b = 1.5pi by definition
 omegab = (varpi_deg - 270.0D0)*radian
 
 ! Eqn 23 of Kipping 2014
 loglikeprior = DLOG( 1.0D0 + eccb*DSIN(omegab) ) - &
                DLOG( 1.0D0 - eccb**2 ) + &
                (beta0 - 1.0D0)*DLOG(1.0D0 - eccb) + &
                (alpha0 - 1.0D0)*DLOG(eccb) - & 
                DLOG(intconst0)

END FUNCTION
! ======================================================================

! ======================================================================
FUNCTION loglikebetaprior(eccc)

 implicit none

 ! Inputs
 REAL(8) :: eccc
 ! Constants
 REAL(8), PARAMETER :: alpha0 = 0.867D0
 REAL(8), PARAMETER :: beta0 = 3.03D0
 REAL(8), PARAMETER :: intconst1 = 0.42718563693158357D0 ! needs to be updated by hand
 ! Output
 REAL(8) :: loglikebetaprior
 
 loglikebetaprior = (beta0 - 1.0D0)*DLOG(1.0D0 - eccc) + &
                    (alpha0 - 1.0D0)*DLOG(eccc) - & 
                    DLOG(intconst1)

END FUNCTION
! ======================================================================

! ======================================================================
FUNCTION chaotic(P1,e1,w1,mu1,P2,e2,w2,mu2)

 ! tests Hadden & Lithwick (2018) chaos criterion
 ! technically 1 is inner planet, 2 is outer planet, but we use an abs
 ! function around the a2-a1 functions which should remove this issue
 ! because cos(x) is symmetric about x too

 implicit none

 ! Inputs
 REAL(8) :: mu1, mu2 ! mp/mstar
 REAL(8) :: e1, e2 ! eccenticities
 REAL(8) :: w1, w2 ! varpi of each planet
 REAL(8) :: P1, P2 ! periods
 ! Intermediate
 REAL(8) :: Zsys, a1, a2, Zcrit
 ! Constants
 REAL(8), PARAMETER :: roothalf = 0.7071067811865476D0
 REAL(8), PARAMETER :: third = 0.33333333D0
 REAL(8), PARAMETER :: twothirds = 0.66666667D0
 REAL(8), PARAMETER :: fourthirds = 1.33333333D0
 ! Output
 INTEGER :: chaotic ! 1 = yes, 0 = no

 ! Z of the system, Eq 16 of Hadden & Lithwick (2018)
 Zsys = roothalf*DSQRT( e1**2 + e2**2 - 2.0D0*e1*e2*DCOS(w1-w2) )

 ! Not strictly true, but OK since we only ever use a1/a2 ratios
 a1 = P1**(twothirds)
 a2 = P2**(twothirds)
 
 ! Critical Z, Eq 19 of Hadden & Lithwick (2018)
 Zcrit = -2.2D0*(mu1+mu2)**third*DABS(a2/(a2-a1))**fourthirds
 Zcrit = roothalf*DABS((a2-a1)/a1)*DEXP( Zcrit )
 
 IF( Zsys .GT. Zcrit ) THEN
   chaotic = 1
 ELSE
   chaotic = 0
 END IF

END FUNCTION
! ======================================================================

! ======================================================================
FUNCTION hillunstable(P1,e1,mu1,P2,e2,mu2)

 ! tests Petit et al. (2018) Hill criterion
 ! planet 1 needs to be the inner

 implicit none

 ! Inputs
 REAL(8) :: mu1, mu2 ! mp/mstar
 REAL(8) :: e1, e2 ! eccenticities
 REAL(8) :: P1, P2 ! periods
 REAL(8) :: inc1, inc2 ! inclinations, unclear how defined TBH
 ! Intermediate
 REAL(8) :: alpha, gamma, epsilon, P1P2third, e1x, e2x
 REAL(8) :: Csys, Ccrit
 ! Constants
 REAL(8), PARAMETER :: k1 = 4.3267487109222245D0 ! 3^(4/3)
 REAL(8), PARAMETER :: third = 0.33333333D0
 REAL(8), PARAMETER :: twothirds = 0.66666667D0
 ! Output
 INTEGER :: hillunstable ! 1 = yes, 0 = no

 gamma = (mu1/mu2)
 P1P2third = (P1/P2)**third
 alpha = P1P2third**2
 epsilon = mu1 + mu2
 e1x = DSQRT( 1.0D0 - e1**2 )
 e2x = DSQRT( 1.0D0 - e2**2 )

 ! C of the system, Eq 23 of Petit et al. (2018)
 inc1 = 0.0D0 ! coplanar approx
 inc2 = 0.0D0 ! coplanar approx
 Csys = gamma*P1P2third*( 1.0D0 - e1x*DCOS(inc1) ) + 1.0D0 - e2x*DCOS(inc2)
 
 ! Critical C, Eq 26 of Petit et al. (2018)
 Ccrit = ( k1*epsilon**twothirds*gamma )/( 1.0D0 + gamma )**2
 Ccrit = (alpha/(gamma+alpha))*(1.0D0+Ccrit)
 Ccrit = gamma*P1P2third + 1.0D0 - DSQRT( (1.0D0+gamma)**3*Ccrit )
 !write(*,*) 'inputs = ',P1,mu1,P2,mu2
 !write(*,*) 'intermediates = ',epsilon,gamma,alpha,P1P2third
 !write(*,*) 'Ccrit = ',Ccrit
 
 IF( Csys .GT. Ccrit ) THEN
   hillunstable = 1
 ELSE
   hillunstable = 0
 END IF

END FUNCTION
! ======================================================================

END MODULE like

