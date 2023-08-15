PROGRAM main

	use params
	use nestwrapper
      
	implicit none
	
	INTEGER :: i, j

      	!no parameters to wrap around
      	nest_pWrap = 0
        nest_pWrap(10) = 1 !phi
        nest_pWrap(11) = 1 !cosis
        nest_pWrap(12) = 1 !Os
        !nest_pWrap(16) = 1 !wp

        ! OBTAIN PRIORS
        !--------------------------------------------------------------------------
        ! Open up jammi_in.txt and read file
        open(1,FILE='example_LRVmoon/jammi_in.jam',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,sdim
          read(1,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
          write(*,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
        END DO
        close(1)
				
        ! OBTAIN BLENDS
        !--------------------------------------------------------------------------
        ! Open up jammi_in.txt and read file
        open(2784,FILE='example_LRVmoon/blends.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,NQ
          read(2784,*) gam(j)
					gam(j) = 1.0D0/gam(j)
        END DO
        close(2784)

        ! Mode 0 override
        DO j=1,sdim
          IF( Rflag(j) .EQ. 0 ) THEN
            Rmin(j) = Rsol(j) !- 1.0D-9
            Rmax(j) = Rsol(j) !+ 1.0D-9
          END IF
        END DO
        
        ! Write out
        DO j=1,sdim
        write(*,*) Rmin(j),Rmax(j)
        END DO
        !--------------------------------------------------------------------------

        ! Read-in Kepler
        IF( priflag .EQ. 1 ) THEN
          write(*,*) 'Reading in primary data...'
          open(unit=30301,file='example_LRVmoon/seriesP.jam')
          DO i = 1,nplen
            read(30301,*) tp(i),fp(i),sigfp(i),epochp_seq(i)
            fpwei(i) = fp(i)/(sigfp(i)**2)
            sigfpwei(i) = 1.0D0/(sigfp(i)**2)
            logsigfp(i) = DLOG(sigfp(i)**2)
          END DO
          close(30301)
        END IF

        ! Read-in seriesR.dat
        IF( rvflag .EQ. 1 ) THEN
          write(*,*) 'Reading in rv data...'
          open(unit=30302,file='example_LRVmoon/seriesR.dat')
          DO i = 1,nrlen
            read(30302,*) tr(i),rv(i),sigrv(i),rvmode(i)
          END DO
          close(30302)
        END IF
				
        ! Read-in TESS
        write(*,*) 'Reading in TESS data...'
        open(unit=30303,file='example_LRVmoon/lightcurves/TESS.dat')
        DO i = 1,nslen
          read(30303,*) ts(i),fs(i),sigfs(i)
					epochs_seq(i) = 9
          fswei(i) = fs(i)/(sigfs(i)**2)
          sigfswei(i) = 1.0D0/(sigfs(i)**2)
          logsigfs(i) = DLOG(sigfs(i)**2)
          ! tsdeviant(i) = time from Rsol ephemeris
          tsdeviant(i) = ts(i) - Rsol(5) &
                         - Rsol(4)*DNINT( (ts(i)-Rsol(5)) / Rsol(4) )
        END DO
        close(30303)
				
        ! Read-in BARON
        write(*,*) 'Reading in BARON data...'
        open(unit=30304,file='example_LRVmoon/lightcurves/BARON.dat')
        DO i = 1,nqlen
          read(30304,*) tq(i),fq(i),sigfq(i)
					epochq_seq(i) = 10
          fqwei(i) = fq(i)/(sigfq(i)**2)
          sigfqwei(i) = 1.0D0/(sigfq(i)**2)
          logsigfq(i) = DLOG(sigfq(i)**2)
          ! tqdeviant(i) = time from Rsol ephemeris
          tqdeviant(i) = tq(i) - Rsol(5) &
                         - Rsol(4)*DNINT( (tq(i)-Rsol(5)) / Rsol(4) )
        END DO
        close(30304)
				
        ! Read-in LCO
        write(*,*) 'Reading in LCO data...'
        open(unit=30305,file='example_LRVmoon/lightcurves/LCO_Teide.dat')
        DO i = 1,ntlen
          read(30305,*) tt(i),ft(i),sigft(i)
					epocht_seq(i) = 11
          ftwei(i) = ft(i)/(sigft(i)**2)
          sigftwei(i) = 1.0D0/(sigft(i)**2)
          logsigft(i) = DLOG(sigft(i)**2)
          ! tqdeviant(i) = time from Rsol ephemeris
          ttdeviant(i) = tt(i) - Rsol(5) &
                         - Rsol(4)*DNINT( (tt(i)-Rsol(5)) / Rsol(4) )
        END DO
        close(30305)
				
        ! Read-in Whitin
        write(*,*) 'Reading in Whitin data...'
        open(unit=30306,file='example_LRVmoon/lightcurves/Whitin.dat')
        DO i = 1,nulen
          read(30306,*) tu(i),fu(i),sigfu(i)
					epochu_seq(i) = 12
          fuwei(i) = fu(i)/(sigfu(i)**2)
          sigfuwei(i) = 1.0D0/(sigfu(i)**2)
          logsigfu(i) = DLOG(sigfu(i)**2)
          ! tqdeviant(i) = time from Rsol ephemeris
          tudeviant(i) = tu(i) - Rsol(5) &
                         - Rsol(4)*DNINT( (tu(i)-Rsol(5)) / Rsol(4) )
        END DO
        close(30306)

      	call nest_Sample
END
