! ----------------------------------------------------------------------
! Copyright: CNRS - Université de Grenoble Alpes
!
! Contributors : Jean-Michel Brankart
!
! Jean-Michel.Brankart@univ-grenoble-alpes.fr
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
!---------------------------------------------------------------------
!
!                        MODULE STORFG
!
!---------------------------------------------------------------------
! Random field generator (Gaussian, with given power spectrum)
! by Jean-Michel Brankart, September 2014
! ----------------------------------------------------------------------

MODULE ensdam_storfg
   !!======================================================================
   !!                       ***  MODULE  storfg  ***
   !! Random field generator (Gaussian, with given power spectrum)
   !!=====================================================================
   !! This module implements the method described in:
   !! Mejia, J. and Rodriguez-Iturbe I. (1974).
   !! On the synthesis of random field sampling from the spectrum:
   !! an application to the generation of hydrological spatial processes.
   !! Water Resources Research. 10(4), 705-711.

   !!----------------------------------------------------------------------
   !! Public routines :
   !!   def_spect_init  : initialize definition of power spectra
   !!   def_spect_power : define power spectra (one by one)
   !!   def_sample_size : define sample size (in 1d and 2d)
   !!   sample_freq_1d  : sample new wave numbers from 1d power spectrum
   !!   sample_freq_2d  : sample new wave numbers from 2d power spectrum
   !!   sample_freq_2s  : sample new wave numbers from 2d power spectrum (on the sphere)
   !!   gen_field_1d    : genrate 1d random field
   !!   gen_field_2d    : genrate 2d random field
   !!   gen_field_2s    : genrate 2d random field (on the sphere)
   USE ensdam_storng
   USE ensdam_sphylm
#if defined MPI
   use mpi
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC def_spect_init
   PUBLIC def_spect_power
   PUBLIC def_sample_size
   PUBLIC sample_freq_1d
   PUBLIC sample_freq_2d
   PUBLIC sample_freq_2s
   PUBLIC gen_field_1d
   PUBLIC gen_field_2d
   PUBLIC gen_field_2s

#if defined MPI
   ! Public definitions for MPI
   INTEGER, PUBLIC, SAVE  :: mpi_comm_storfg=mpi_comm_world   ! definition of module global communicator
   INTEGER, save :: mpi_code
#endif

   ! Resolution for precomputation of Legendre polynomials (in degrees)
   REAL(KIND=8), PUBLIC, SAVE :: storfg_ylm_resolution=0.01

   INTEGER :: jpfreq    ! number of wave numbers defining power spectra
   INTEGER :: jpspct1d  ! number of 1d power spectra
   INTEGER :: jpspct2d  ! number of 2d power spectra
   INTEGER :: jpspct2s  ! number of 2d power spectra, on the sphere
   INTEGER :: jpsmp1d   ! size of frequency sample (in 1d)
   INTEGER :: jpsmp2d   ! size of frequency sample (in 2d)
   INTEGER :: jpsmp2s   ! size of frequency sample (in 2d, on the sphere)

   REAL(KIND=8), PARAMETER :: twopi=2*3.1415926535897932384626
   REAL(KIND=8), PARAMETER :: deg2rad=twopi/360.

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spect_freq_1d  ! list of wave numbers (1d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spect_freq_2d  ! list of wave numbers (2d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spect_pow_1d   ! 1d power spectra
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spect_pow_2d   ! 2d power spectra
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: spect_pow_2s   ! 2d power spectra, on the sphere

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_freq_1d    ! wave number samples (1d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_phas_1d    ! phase samples (1d)

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_freq_2d    ! wave number samples (2d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_phas_2d    ! phase samples (2d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_dir1_2d    ! dir (comp 1) samples (2d)
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: smp_dir2_2d    ! dir (comp 2) samples (2d)

   INTEGER,      DIMENSION(:,:),   ALLOCATABLE :: smp_freq_2s  ! wave number samples (2s)
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: smp_dir_2s   ! direction samples (2s)
   LOGiCAL :: pleg_precomputed = .FALSE.

   INTERFACE gen_field_2s
     MODULE PROCEDURE gen_field_2s_1d, gen_field_2s_2d, gen_field_2s_new, gen_field_2s_new2
   END INTERFACE

   INTERFACE sample_freq_2s
     MODULE PROCEDURE sample_freq_2s_tabulated, sample_freq_2s_function
   END INTERFACE

   INTERFACE
     FUNCTION fun_spect_pow_2s(l)      ! callback function for power spectrum on the sphere
     IMPLICIT NONE
     INTEGER, intent(in) :: l
     REAL(KIND=8) :: fun_spect_pow_2s
     END FUNCTION
   END INTERFACE

   INTERFACE
     FUNCTION fun_pow_spect_sph(l,m)   ! callback function for power spectrum on the sphere
     IMPLICIT NONE
     INTEGER, intent(in) :: l,m
     REAL(KIND=8) :: fun_pow_spect_sph
     END FUNCTION
   END INTERFACE

   ! Public definitions needed by python/julia APIs
   PUBLIC gen_field_2s_new, gen_field_2s_new2, fun_pow_spect_sph

CONTAINS

   SUBROUTINE def_spect_init( kjpfreq, kjpspct1d, kjpspct2d, kjpspct2s )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE def_spect_init  ***
      !! 
      !! ** Purpose :   define number of spectra and wave numbers
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kjpfreq   ! number of wave numbers
      INTEGER, INTENT( in ) ::   kjpspct1d ! number of 1d spectra
      INTEGER, INTENT( in ) ::   kjpspct2d ! number of 2d spectra
      INTEGER, INTENT( in ) ::   kjpspct2s ! number of 2d spectra, on the sphere

      jpfreq   = kjpfreq
      jpspct1d = kjpspct1d
      jpspct2d = kjpspct2d
      jpspct2s = kjpspct2s

      ! Allocate list of wave numbers
      IF (jpspct1d>0) ALLOCATE(spect_freq_1d(jpfreq,jpspct1d))
      IF (jpspct2d>0) ALLOCATE(spect_freq_2d(jpfreq,jpspct2d))

      ! Allocate power spectra
      IF (jpspct1d>0) ALLOCATE(spect_pow_1d(jpfreq,jpspct1d))
      IF (jpspct2d>0) ALLOCATE(spect_pow_2d(jpfreq,jpspct2d))
      IF (jpspct2s>0) ALLOCATE(spect_pow_2s(jpfreq,jpspct2s))

   END SUBROUTINE def_spect_init


   SUBROUTINE def_spect_power(ktyp,kspct,kspect_freq,kspect_power)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE def_spect_power  ***
      !! 
      !! ** Purpose :   define power spectra (one by one)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  ktyp   ! type of power sepctrum to define (1d, 2d or 2s)
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kspect_freq  ! wave numbers
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kspect_power ! power spectrum

      ! Define spectra for 1d fields
      IF (ktyp==1) THEN
        IF ( (kspct<1) .OR. (kspct>jpspct1d) ) STOP 'Bad sepctrum index in def_spect_power'
        spect_freq_1d(:,kspct) = kspect_freq(:)
        spect_pow_1d(:,kspct) = kspect_power(:)
      ENDIF

      ! Define spectra for 2d fields
      IF (ktyp==2) THEN
        IF ( (kspct<1) .OR. (kspct>jpspct2d) ) STOP 'Bad sepctrum index in def_spect_power'
        spect_freq_2d(:,kspct) = kspect_freq(:)
        spect_pow_2d(:,kspct) = kspect_power(:)
      ENDIF

      ! Define spectra for 2d fields
      IF (ktyp==3) THEN
        IF ( (kspct<1) .OR. (kspct>jpspct2s) ) STOP 'Bad sepctrum index in def_spect_power'
        spect_pow_2s(:,kspct) = kspect_power(:)
      ENDIF

   END SUBROUTINE def_spect_power


   SUBROUTINE def_sample_size( kjpsmp1d, kjpsmp2d, kjpsmp2s, kseed )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE def_sample_size
      !! 
      !! ** Purpose :   define size of frequency samples
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: kjpsmp1d ! size of frequancey sample for 1d spectra
      INTEGER, INTENT( in ) :: kjpsmp2d ! size of frequancey sample for 2d spectra
      INTEGER, INTENT( in ) :: kjpsmp2s ! size of frequancey sample for 2s spectra
      INTEGER, INTENT( in ) :: kseed    ! index of seed for RNG
      !!
      INTEGER :: jseed
      INTEGER(KIND=8) :: zseed1, zseed2, zseed3, zseed4

      jpsmp1d = kjpsmp1d
      jpsmp2d = kjpsmp2d
      jpsmp2s = kjpsmp2s

      ! Allocate samples for 1d fields
      IF (jpspct1d>0) THEN
        ALLOCATE(smp_freq_1d(jpsmp1d,jpspct1d))
        ALLOCATE(smp_phas_1d(jpsmp1d,jpspct1d))
      ENDIF

      ! Allocate samples for 2d fields
      IF (jpspct2d>0) THEN
        ALLOCATE(smp_freq_2d(jpsmp2d,jpspct2d))
        ALLOCATE(smp_phas_2d(jpsmp2d,jpspct2d))
        ALLOCATE(smp_dir1_2d(jpsmp2d,jpspct2d))
        ALLOCATE(smp_dir2_2d(jpsmp2d,jpspct2d))
      ENDIF

      ! Allocate samples for 2d fields
      IF (jpspct2s>0) THEN
        ALLOCATE(smp_freq_2s(jpsmp2s,jpspct2s))
        ALLOCATE(smp_dir_2s(jpsmp2s,-jpfreq:jpfreq,jpspct2s))
      ENDIF

      ! Seed random number generator
      CALL kiss_reset( )
      DO jseed = 1, kseed
        zseed1 = kiss() ; zseed2 = kiss() ; zseed3 = kiss() ; zseed4 = kiss()
      END DO
      CALL kiss_seed( zseed1, zseed2, zseed3, zseed4 )

   END SUBROUTINE def_sample_size


   SUBROUTINE sample_freq_1d(kspct)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sample_freq_1d  ***
      !! 
      !! ** Purpose :   sample new set of wave numbers from 1d power spectrum
      !!                sample corresponding random phase, uniformly in [0,2pi]
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      !!
      INTEGER :: jsmp, jfreq
      LOGICAL :: rejected
      REAL(KIND=8) :: f, fmin, fmax, pow, pow_min, pow_max, pow_avg
      REAL(KIND=8) :: uran   ! Uniform random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: eps = 0.001  ! Small number

      ! Get maximum and minimum power spectrum
      pow_max = MAXVAL(spect_pow_1d(:,kspct))
      pow_min = eps * pow_max

      ! Get maximum wave number
      fmax = spect_freq_1d(1,kspct)
      DO jfreq = 1,jpfreq
        IF (spect_pow_1d(jfreq,kspct).GT.pow_min) THEN
          fmax = spect_freq_1d(jfreq,kspct)
        ENDIF
      ENDDO

      ! Get minimum wave number
      fmin = spect_freq_1d(jpfreq,kspct)
      DO jfreq = jpfreq,1,-1
        IF (spect_pow_1d(jfreq,kspct).GT.pow_min) THEN
          fmin = spect_freq_1d(jfreq,kspct)
        ENDIF
      ENDDO

      ! Sample wave numbers from 1d power spectrum
      DO jsmp = 1, jpsmp1d
        rejected=.TRUE.
        DO WHILE(rejected)
          ! sample uniformly in [fmin,fmax]
          CALL kiss_uniform(uran)
          f = fmin + (fmax-fmin) * uran
          ! interpolate power spectrum at f
          CALL pow_interp(pow,f,spect_freq_1d(:,kspct),spect_pow_1d(:,kspct))
          ! accept with probbaility pow/pow_max
          CALL kiss_uniform(uran)
          rejected = (uran * pow_max) .GT. pow
        ENDDO
        smp_freq_1d(jsmp,kspct) = f
      ENDDO

      ! Sample random phase, uniformly in [0,2pi]
      DO jsmp = 1, jpsmp1d
        CALL kiss_uniform(uran)
        smp_phas_1d(jsmp,kspct) = uran * twopi
      ENDDO

   END SUBROUTINE sample_freq_1d


   SUBROUTINE sample_freq_2d(kspct)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sample_freq_2d  ***
      !! 
      !! ** Purpose :   sample new set of wave numbers from 2d power spectrum
      !!                sample corresponding random phase, uniformly in [0,2pi]
      !!                sample driection, uniformly in [0,2pi]
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      !!
      INTEGER :: jsmp, jfreq
      LOGICAL :: rejected
      REAL(KIND=8) :: f, fmin, fmax, pow, pow_min, pow_max, pow_avg
      REAL(KIND=8) :: uran   ! Uniform random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: eps = 0.001  ! Small number

      ! Get maximum and minimum power spectrum
      pow_max = MAXVAL(spect_pow_2d(:,kspct))
      pow_min = eps * pow_max

      ! Get maximum wave number
      fmax = spect_freq_2d(1,kspct)
      DO jfreq = 1,jpfreq
        IF (spect_pow_2d(jfreq,kspct).GT.pow_min) THEN
          fmax = spect_freq_2d(jfreq,kspct)
        ENDIF
      ENDDO

      ! Get minimum wave number
      fmin = spect_freq_2d(jpfreq,kspct)
      DO jfreq = jpfreq,1,-1
        IF (spect_pow_2d(jfreq,kspct).GT.pow_min) THEN
          fmin = spect_freq_2d(jfreq,kspct)
        ENDIF
      ENDDO

      ! Sample wave numbers from 2d power spectrum
      DO jsmp = 1, jpsmp2d
        rejected=.TRUE.
        DO WHILE(rejected)
          ! sample uniformly in [fmin,fmax]
          CALL kiss_uniform(uran)
          f = fmin + (fmax-fmin) * uran
          ! interpolate power spectrum at f
          CALL pow_interp(pow,f,spect_freq_2d(:,kspct),spect_pow_2d(:,kspct))
          ! accept with probbaility pow/pow_max
          CALL kiss_uniform(uran)
          rejected = (uran * pow_max) .GT. pow
        ENDDO
        smp_freq_2d(jsmp,kspct) = f
      ENDDO

      ! Sample random phase, uniformly in [0,2pi]
      DO jsmp = 1, jpsmp2d
        CALL kiss_uniform(uran)
        smp_phas_2d(jsmp,kspct) = uran * twopi
      ENDDO

      ! Sample direction, uniformly in [0,2pi]
      DO jsmp = 1, jpsmp2d
        CALL kiss_uniform(uran)
        smp_dir1_2d(jsmp,kspct) = COS( uran * twopi )
        smp_dir2_2d(jsmp,kspct) = SIN( uran * twopi )
      ENDDO

   END SUBROUTINE sample_freq_2d


   SUBROUTINE sample_freq_2s_tabulated(kspct)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sample_freq_2s  ***
      !! 
      !! ** Purpose :   sample new set of wave numbers
      !!                   from 2d power spectrum (on the sphere)
      !!                sample corresponding uniformly on the sphere
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      !!
      INTEGER :: jsmp, jfreq, l, lmin, lmax, m
      LOGICAL :: rejected
      REAL(KIND=8) :: pow, pow_min, pow_max, pow_avg, norm
      REAL(KIND=8) :: uran   ! Uniform random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: gran   ! Gaussian random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: eps = 0.001  ! Small number
      REAL(KIND=8) :: latmin, latmax, dlatmax

      ! Get maximum and minimum power spectrum
      pow_max = MAXVAL(spect_pow_2s(:,kspct))
      pow_min = eps * pow_max

      ! Get maximum wave number
      lmax = 0
      DO jfreq = 0,jpfreq-1
        IF (spect_pow_2s(jfreq+1,kspct).GT.pow_min) THEN
          lmax = jfreq
        ENDIF
      ENDDO

      ! Get minimum wave number
      lmin = jpfreq-1
      DO jfreq = jpfreq-1,0,-1
        IF (spect_pow_2s(jfreq+1,kspct).GT.pow_min) THEN
          lmin = jfreq
        ENDIF
      ENDDO

      ! Initialize computation of spherical harmonics
      latmin = -90. ; latmax = 90. ; dlatmax = 0.1
      CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )

      ! Sample wave numbers from power spectrum
      DO jsmp = 1, jpsmp2s
        rejected=.TRUE.
        DO WHILE(rejected)
          ! sample uniformly in [lmin,lmax]
          CALL kiss_uniform(uran)
          l = INT( lmin + (lmax-lmin+1) * uran )
          ! get power spectrum at l
          pow=spect_pow_2s(l+1,kspct)
          ! accept with probbaility pow/pow_max
          CALL kiss_uniform(uran)
          rejected = (uran * pow_max) .GT. pow
        ENDDO
        smp_freq_2s(jsmp,kspct) = l
      ENDDO

      ! Sample direction, uniformly on a sphere of dimension 2l+1
      DO jsmp = 1, jpsmp2d
        ! Sample Gaussian vector in N(0,I), in space of direction 2l+1
        smp_dir_2s(jsmp,:,kspct) = 0.
        DO m = -smp_freq_2s(jsmp,kspct), smp_freq_2s(jsmp,kspct)
          CALL kiss_gaussian(gran)
          smp_dir_2s(jsmp,m,kspct) = gran
        ENDDO
        ! Compute l2-norm of this random vector
        norm = 0.
        DO m = -smp_freq_2s(jsmp,kspct), smp_freq_2s(jsmp,kspct)
          norm = norm + smp_dir_2s(jsmp,m,kspct) * smp_dir_2s(jsmp,m,kspct) 
        ENDDO
        ! Normalize random vector to get it on the unit sphere
        smp_dir_2s(jsmp,:,kspct) = smp_dir_2s(jsmp,:,kspct) / SQRT(norm)
      ENDDO

   END SUBROUTINE sample_freq_2s_tabulated

   SUBROUTINE sample_freq_2s_function(fun_spect,lmin,lmax,kspct)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sample_freq_2s  ***
      !! 
      !! ** Purpose :   sample new set of wave numbers
      !!                   from 2d power spectrum (on the sphere)
      !!                sample corresponding direction uniformly on the sphere
      !!----------------------------------------------------------------------
      PROCEDURE(fun_spect_pow_2s) :: fun_spect ! power spectrum function
      INTEGER, INTENT( in ) :: lmin  ! minimum degree of the spherical harmonics
      INTEGER, INTENT( in ) :: lmax  ! maximum degree of the spherical harmonics
      INTEGER, INTENT( in ) :: kspct ! index of table used to store random draws
      !!
      INTEGER :: jsmp, l, m
      LOGICAL :: rejected
      REAL(KIND=8) :: pow, pow_max, norm
      REAL(KIND=8) :: uran   ! Uniform random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: gran   ! Gaussian random number (forced KIND=8 as in kiss)
      REAL(KIND=8) :: latmin, latmax, dlatmax

      ! Get maximum of power spectrum
      pow_max = 0.0
      DO l = lmin,lmax
        pow = fun_spect(l)
        IF (pow.GT.pow_max) pow_max = pow
      ENDDO

      ! Initialize computation of spherical harmonics
      latmin = -90. ; latmax = 90. ; dlatmax = 0.1
      CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )

      ! Sample wave numbers from power spectrum
      DO jsmp = 1, jpsmp2s
        rejected=.TRUE.
        DO WHILE(rejected)
          ! sample uniformly in [lmin,lmax]
          CALL kiss_uniform(uran)
          l = INT( lmin + (lmax-lmin+1) * uran )
          ! get power spectrum at l
          pow=fun_spect(l)
          ! accept with probbaility pow/pow_max
          CALL kiss_uniform(uran)
          rejected = (uran * pow_max) .GT. pow
        ENDDO
        smp_freq_2s(jsmp,kspct) = l
      ENDDO

      ! Sample direction, uniformly on a sphere of dimension 2l+1
      DO jsmp = 1, jpsmp2s
        ! Sample Gaussian vector in N(0,I), in space of direction 2l+1
        smp_dir_2s(jsmp,:,kspct) = 0.
        DO m = -smp_freq_2s(jsmp,kspct), smp_freq_2s(jsmp,kspct)
          CALL kiss_gaussian(gran)
          smp_dir_2s(jsmp,m,kspct) = gran
        ENDDO
        ! Compute l2-norm of this random vector
        norm = 0.
        DO m = -smp_freq_2s(jsmp,kspct), smp_freq_2s(jsmp,kspct)
          norm = norm + smp_dir_2s(jsmp,m,kspct) * smp_dir_2s(jsmp,m,kspct)
        ENDDO
        ! Normalize random vector to get it on the unit sphere
        smp_dir_2s(jsmp,:,kspct) = smp_dir_2s(jsmp,:,kspct) / SQRT(norm)
      ENDDO

#if defined MPI
      CALL mpi_bcast(smp_freq_2s(:,kspct),jpsmp2s,mpi_double_precision,0,mpi_comm_storfg,mpi_code)
      CALL mpi_bcast(smp_dir_2s(:,:,kspct),jpsmp2s*(2*jpfreq+1),mpi_double_precision,0,mpi_comm_storfg,mpi_code)
#endif

   END SUBROUTINE sample_freq_2s_function

   SUBROUTINE gen_field_1d(kspct,ranfield,x)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_1d  ***
      !! 
      !! ** Purpose :   Generate 1d random field with the required spectrum
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      REAL(KIND=8), DIMENSION(:), INTENT( in )  :: x        ! spatial coordinates
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ranfield ! random field
      !!
      INTEGER :: jsmp

      ranfield(:) = 0.0
      DO jsmp = 1, jpsmp1d
         ranfield(:) = ranfield(:) +                               &
                             COS( smp_freq_1d(jsmp,kspct) * x(:)   &
                                  + smp_phas_1d(jsmp,kspct) )
      ENDDO
      ranfield(:) = ranfield(:) * SQRT(2.0/REAL(jpsmp1d,8))

   END SUBROUTINE gen_field_1d


   SUBROUTINE gen_field_2d(kspct,ranfield,x,y)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_2d  ***
      !! 
      !! ** Purpose :   Generate 2d random field with the required spectrum
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      REAL(KIND=8), DIMENSION(:,:), INTENT( in )  :: x,y      ! spatial coordinates
      REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: ranfield ! random field
      !!
      INTEGER :: jsmp

      ranfield(:,:) = 0.0
      DO jsmp = 1, jpsmp2d
         ranfield(:,:) = ranfield(:,:) +                                   &
                             COS( smp_freq_2d(jsmp,kspct)                  &
                                  * ( smp_dir1_2d(jsmp,kspct) * x(:,:) +   &
                                      smp_dir2_2d(jsmp,kspct) * y(:,:) )   &
                                  + smp_phas_2d(jsmp,kspct) )
      ENDDO
      ranfield(:,:) = ranfield(:,:) * SQRT(2.0/REAL(jpsmp2d,8))

  END SUBROUTINE gen_field_2d

  SUBROUTINE gen_field_2s_2d(kspct,ranfield,lon,lat)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_2s_2d  ***
      !! 
      !! ** Purpose :   Generate 2d random field on the sphere
      !!                    with the required spectrum
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      REAL(KIND=8), DIMENSION(:,:), INTENT( in )  :: lon, lat ! spatial coordinates
      REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: ranfield ! random field
      !!
      INTEGER :: jsmp, ji, jj, l, m
      REAL(KIND=8) :: rcomp

      ! Construct random field by addinp up spherical harmonics
      DO jj = 1, SIZE(ranfield,2)
        CALL gen_field_2s_1d(kspct,ranfield(:,jj),lon(:,jj),lat(:,jj))
      ENDDO

  END SUBROUTINE gen_field_2s_2d

  SUBROUTINE gen_field_2s_1d(kspct,ranfield,lon,lat)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_2s_1d  ***
      !! 
      !! ** Purpose :   Generate random field on the sphere
      !!                   with the required spectrum
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kspct  ! index of power spectrum to define
      REAL(KIND=8), DIMENSION(:), INTENT( in )  :: lon, lat ! spatial coordinates
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ranfield ! random field
      !!
      INTEGER :: jsmp, ji, l, m
      REAL(KIND=8) :: rcomp

      ! Construct random field by addinp up spherical harmonics
      ranfield(:) = 0.0
      DO ji = 1, SIZE(ranfield,1)
        ! Loop on (random) degree of spherical harmonics
        DO jsmp = 1, jpsmp2s
          l = smp_freq_2s(jsmp,kspct)
          ! Loop on order of spherical harmonics
          DO m = -l, l
            rcomp = smp_dir_2s(jsmp,m,kspct)
            rcomp = rcomp * ylm(l,m,lon(ji),lat(ji))
            ranfield(ji) = ranfield(ji) + rcomp
          ENDDO
        ENDDO
      ENDDO
      ranfield(:) = ranfield(:) / SQRT(REAL(jpsmp2s,8))

  END SUBROUTINE gen_field_2s_1d

  SUBROUTINE gen_field_2s_new(ranfield,lon,lat,pow_spect,lmin,lmax)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_2s_new  ***
      !! 
      !! ** Purpose :   Generate random field on the sphere
      !!                   with the required spectrum
      !!
      !!----------------------------------------------------------------------
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ranfield ! random field
      REAL(KIND=8), DIMENSION(:), INTENT( in )  :: lon, lat ! spatial coordinates
      PROCEDURE(fun_pow_spect_sph) :: pow_spect ! square root of power spectrum
      INTEGER, INTENT( in ) ::  lmin, lmax  ! range of degrees
      !!
      INTEGER :: jsmp, ji, jpi, l, m
      REAL(KIND=8) :: gran, amplitude
      REAL(KIND=8) :: latmin, latmax, dlatmax

      jpi = SIZE(ranfield,1)

      ! Initialize computation of spherical harmonics
      IF (.NOT.pleg_precomputed) THEN
        latmin = -90 ; latmax = 90 ; dlatmax = storfg_ylm_resolution
        CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )
        pleg_precomputed = .TRUE.
      ENDIF

      ! Construct random field by addinp up spherical harmonics
      ranfield(:) = 0.0
      ! Loop on (random) degree of spherical harmonics
      DO l = lmin, lmax
        ! Loop on order of spherical harmonics
        DO m = -l, l
          amplitude = SQRT(pow_spect(l,m))
          CALL kiss_gaussian(gran)
#if defined MPI
          CALL mpi_bcast(gran,1,mpi_double_precision,0,mpi_comm_storfg,mpi_code)
#endif
          DO ji = 1, jpi
            ranfield(ji) = ranfield(ji) + gran * amplitude * ylm(l,m,lon(ji),lat(ji))
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE gen_field_2s_new

  SUBROUTINE gen_field_2s_new2(ranfield,lon,lat,pow_spect,lmin,lmax)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE gen_field_2s_new  ***
      !! 
      !! ** Purpose :   Generate random field on the sphere
      !!                   with the required spectrum
      !!
      !!----------------------------------------------------------------------
      REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: ranfield ! random field
      REAL(KIND=8), DIMENSION(:,:), INTENT( in )  :: lon, lat ! spatial coordinates
      PROCEDURE(fun_pow_spect_sph) :: pow_spect ! square root of power spectrum
      INTEGER, INTENT( in ) ::  lmin, lmax  ! range of degrees
      !!
      INTEGER :: jsmp, ji, jj, jpi, jpj, l, m
      REAL(KIND=8) :: gran, amplitude
      REAL(KIND=8) :: latmin, latmax, dlatmax

      jpi = SIZE(ranfield,1)
      jpj = SIZE(ranfield,2)

      ! Initialize computation of spherical harmonics
      IF (.NOT.pleg_precomputed) THEN
        latmin = -90 ; latmax = 90 ; dlatmax = storfg_ylm_resolution
        CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )
        pleg_precomputed = .TRUE.
      ENDIF

      ! Construct random field by addinp up spherical harmonics
      ranfield(:,:) = 0.0
      ! Loop on (random) degree of spherical harmonics
      DO l = lmin, lmax
        ! Loop on order of spherical harmonics
        DO m = -l, l
          amplitude = SQRT(pow_spect(l,m))
          CALL kiss_gaussian(gran)
#if defined MPI
          CALL mpi_bcast(gran,1,mpi_double_precision,0,mpi_comm_storfg,mpi_code)
#endif
          DO ji = 1, jpi
          DO jj = 1, jpj
            ranfield(ji,jj) = ranfield(ji,jj) + gran * amplitude * ylm(l,m,lon(ji,jj),lat(ji,jj))
          ENDDO
          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE gen_field_2s_new2

  SUBROUTINE pow_interp(pow,f,ftab,powtab)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE pow_interp  ***
      !! 
      !! ** Purpose :   Interpolate in power scpetrum at a given wave number
      !!
      !!----------------------------------------------------------------------
      REAL(KIND=8), INTENT( in )  :: f       ! input wave number
      REAL(KIND=8), INTENT( out ) :: pow     ! output power value
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ftab,powtab ! power spectrum
      !!
      INTEGER :: n, i, il, ih, im
      REAL(KIND=8) :: r

      ! Size of wave number array
      n = SIZE(ftab,1)

      ! Dichotomy in wave number array
      il = 1
      ih = n
      DO WHILE (il.LT.ih-1)
        im = (il+ih)/2
        IF (f.GT.ftab(im)) THEN
          il = im
        ELSE
          ih = im
        ENDIF
      ENDDO

      ! Linear interpolation
      r = ( f - ftab(il) ) / ( ftab(il+1) - ftab(il) )
      pow = powtab(il) + ( powtab(il+1)-powtab(il) ) * r

  END SUBROUTINE pow_interp

END MODULE ensdam_storfg
