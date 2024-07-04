!---------------------------------------------------------------------
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
!                        MODULE SPHYLM
!
!---------------------------------------------------------------------
! Projection of ocean fields on spherical harmonics
! by Jean-Michel Brankart, December 2017
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! init_ylm : initialize computation of spherical harmonics
! init_regr_ylm : initialize regression of observations
! proj_ylm : project on spherical harmonics
! back_ylm : transform back on the sphere
! regr_ylm : regression of observations on spherical harmonics
! disp_ylm : output one single spherical harmonics
! ylm : evaluate spherical harmonics
! ----------------------------------------------------------------------
MODULE ensdam_sphylm
#if defined MPI
        use mpi
#endif
        implicit none
        private
#ifdef MPI_INCLUDE
      include "mpif.h"
#endif

        public init_ylm
        public init_regr_ylm
        public proj_ylm
        public back_ylm
        public back_ylm_loc
        public regr_ylm
        public disp_ylm
        public ylm

        ! Private variables/parameters 
        INTEGER, save :: jpl     ! maximum degree of Legendre polynomials
        INTEGER, save :: jplat  ! number latitude points in saved Legendre polynomials
        INTEGER, save :: jlmin  ! minimum degree of Legendre polynomials

        REAL(KIND=8), save :: lat0, dlat
        REAL(KIND=8), dimension(:,:,:), allocatable, save :: pleg   ! Legendre polynomials

        REAL(KIND=8), parameter :: twopi=2*3.1415926535897932384626
        REAL(KIND=8), parameter :: deg2rad=twopi/360.

        ! Public variables ruling the regression of observations
        CHARACTER(len=80), public, save :: regr_type='local' ! regression type
        INTEGER, public, save :: regr_maxiter=50   ! maximum number of iteration
        INTEGER, public, save :: regr_maxbloc=1    ! maximum size of local blocks
        INTEGER, public, save :: regr_overlap=1    ! overlapping of local blocks
        REAL(KIND=8), public, save :: regr_epsilon=0.01 ! relative difference for convergence test
        REAL(KIND=8), public, save :: regr_rho=1.0 ! weight according to signal std (rho=1) or not (rho=0)

        ! Variables for parallel computation
        LOGICAL, public, save :: external_vector_decomposition=.FALSE. ! external vector decomposition
        INTEGER, save :: jpproc=1   ! number of processors
        INTEGER, save :: jproc=0    ! index of current processor
#if defined MPI
        INTEGER, public, save :: mpi_comm_sphylm=mpi_comm_world   ! definition of module global communicator
        INTEGER, save :: mpi_code
#endif

      CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE init_ylm( kjpl, kjlmin, latmin, latmax, dlatmax )
!----------------------------------------------------------------------
!                  ***  init_ylm  ***
! 
! ** Purpose :   initialize computation of spherical harmonics
! ** Method  :   precompute Legendre polynomials
!----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT( in ) :: kjpl, kjlmin
        REAL(KIND=8), INTENT( in ) :: latmin, latmax, dlatmax

        INTEGER :: allocok, jlat
        REAL(KIND=8) :: xlat, lat
        LOGICAL :: correctshape

        jpl = kjpl
        jlmin = kjlmin
        jplat = INT( (latmax - latmin) / dlatmax )

        lat0 = latmin
        dlat = (latmax - latmin) / jplat

! Allocate Legendre polynomials array
        IF (allocated(pleg)) THEN
          correctshape = (size(pleg,1).EQ.(jplat+1)) 
          correctshape = correctshape .AND. (size(pleg,2).EQ.(jpl+2)) 
          correctshape = correctshape .AND. (size(pleg,3).EQ.(jpl+2)) 
          IF (.NOT.correctshape) deallocate(pleg)
        ENDIF

        IF (.NOT.allocated(pleg)) THEN
          allocate(pleg(0:jplat,0:jpl+1,0:jpl+1), stat=allocok)
          IF (allocok.NE.0) STOP 'Allocation error in init_ylm'
        ENDIF

! Compute Legendre polynomials
        DO jlat = 0,jplat
          lat = latmin + jlat * dlat
          IF (lat>latmax) lat = latmax
          xlat = COS( (90.-lat) * deg2rad )
          CALL plm(xlat,jlat)
        ENDDO

#if defined MPI
        IF (external_vector_decomposition) THEN
          jpproc=1 ; jproc=0
        ELSE
          CALL mpi_comm_size(mpi_comm_sphylm,jpproc,mpi_code)
          CALL mpi_comm_rank(mpi_comm_sphylm,jproc,mpi_code)
        ENDIF
#endif

        END SUBROUTINE init_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE init_regr_ylm( ktype, kmaxiter, kmaxbloc, koverlap, kepsilon, krho)
!----------------------------------------------------------------------
!                  ***  init_regr_ylm  ***
! 
! ** Purpose :   set parameters for regression on spherical harmonics
!----------------------------------------------------------------------
        IMPLICIT NONE
        CHARACTER(len=*), INTENT( in ) :: ktype
        INTEGER, INTENT( in ) :: kmaxiter, kmaxbloc, koverlap
        REAL(KIND=8), INTENT( in ) :: kepsilon, krho

        IF (kmaxiter.LT.1) THEN
          PRINT *, 'Error in init_regr_ylm: maxiter=',kmaxiter
          STOP 'Stopping now'
        ENDIF

        IF (kmaxbloc.LT.1) THEN
          PRINT *, 'Error in init_regr_ylm: maxbloc=',kmaxbloc
          STOP 'Stopping now'
        ENDIF

        IF (koverlap.LT.1) THEN
          PRINT *, 'Error in init_regr_ylm: overlap=',koverlap
          STOP 'Stopping now'
        ENDIF

        IF ((kepsilon.LE.0.0).OR.(kepsilon.GE.1.0)) THEN
          PRINT *, 'Error in init_regr_ylm: epsilon=',kepsilon
          STOP 'Stopping now'
        ENDIF

        IF ((krho.LT.0.0).OR.(krho.GT.1.0)) THEN
          PRINT *, 'Error in init_regr_ylm: rho=',krho
          STOP 'Stopping now'
        ENDIF

        regr_type=ktype
        regr_maxiter=kmaxiter
        regr_maxbloc=kmaxbloc
        regr_overlap=koverlap
        regr_epsilon=kepsilon
        regr_rho=krho

        RETURN

        END SUBROUTINE init_regr_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE proj_ylm( kproj, ktab, klon, klat )
!----------------------------------------------------------------------
!                  ***  proj_ylm  ***
! 
! ** Purpose :   project on spherical harmonics of degree l and order m
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ktab, klon, klat
        REAL(KIND=8), DIMENSION(0:,-jpl:), INTENT( out ) :: kproj

        INTEGER :: jpi, ji, k, l, m

        jpi=SIZE(ktab,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in proj_ylm'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in proj_ylm'
        IF (UBOUND(kproj,1).NE.jpl) STOP 'Inconsistent size in proj_ylm'
        IF (UBOUND(kproj,2).NE.jpl) STOP 'Inconsistent size in proj_ylm'

        ! Initialize to compute the whole integral inside the routine
        kproj(:,:)=0.0

        ! Loop on spherical harmonics
        ! (reduced to one single loop for an easy parallelization)
        DO k=jproc+jlmin*jlmin,(jpl+1)*(jpl+1)-1,jpproc
          l=INT(SQRT(REAL(k,8)))
          m=k-l*l-l
          
          DO ji=1,jpi
            kproj(l,m) = kproj(l,m) + ktab(ji)                     &  
     &                     * ylm(l,m,klon(ji),klat(ji))
          ENDDO
        ENDDO

#if defined MPI
        IF (.NOT.external_vector_decomposition) THEN
          CALL MPI_ALLREDUCE (MPI_IN_PLACE, kproj, (jpl+1)*(2*jpl+1),  &
     &                        MPI_DOUBLE_PRECISION,                    &
     &                        MPI_SUM,mpi_comm_sphylm,mpi_code)
        ENDIF
#endif

        END SUBROUTINE proj_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE back_ylm( kproj, ktab, klon, klat )
!----------------------------------------------------------------------
!                  ***  init_ylm  ***
! 
! ** Purpose :   transform back on the sphere
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(0:,-jpl:), INTENT( in ) :: kproj
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ktab

        INTEGER :: jpi, ji, k, l, m

        jpi=SIZE(ktab,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in back_ylm'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in back_ylm'
        IF (UBOUND(kproj,1).NE.jpl) STOP 'Inconsistent size in back_ylm'
        IF (UBOUND(kproj,2).NE.jpl) STOP 'Inconsistent size in back_ylm'

        ktab(:)=0.0
        ! Loop on spherical harmonics
        ! (reduced to one single loop for an easy parallelization)
        DO k=jproc+jlmin*jlmin,(jpl+1)*(jpl+1)-1,jpproc
          l=INT(SQRT(REAL(k,8)))
          m=k-l*l-l
          
          DO ji=1,jpi
            ktab(ji) = ktab(ji) + kproj(l,m)                   &
     &                     * ylm(l,m,klon(ji),klat(ji))
          ENDDO
        ENDDO

#if defined MPI
        IF (.NOT.external_vector_decomposition) THEN
          CALL MPI_ALLREDUCE (MPI_IN_PLACE, ktab, jpi,           &
     &                        MPI_DOUBLE_PRECISION,              &
     &                        MPI_SUM,mpi_comm_sphylm,mpi_code)
        ENDIF
#endif

        END SUBROUTINE back_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE back_ylm_loc( kproj, ktab, klon, klat, kl0, kl1 )
!----------------------------------------------------------------------
!                  ***  back_ylm_loc  ***
! 
! ** Purpose :   transform back on the sphere
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(0:,-jpl:), INTENT( in ) :: kproj
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ktab
        INTEGER, INTENT( in ) :: kl0, kl1

        INTEGER :: jpi, ji, k, l, m

        jpi=SIZE(ktab,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in back_ylm_loc'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in back_ylm_loc'
        IF (UBOUND(kproj,1).NE.jpl) STOP 'Inconsistent size in back_ylm_loc'
        IF (UBOUND(kproj,2).NE.jpl) STOP 'Inconsistent size in back_ylm_loc'

        ktab(:)=0.0
        ! Loop on spherical harmonics
        ! (reduced to one single loop for an easy parallelization)
        DO k=jproc+kl0*kl0,(kl1+1)*(kl1+1)-1,jpproc
          l=INT(SQRT(REAL(k,8)))
          m=k-l*l-l
          
          DO ji=1,jpi
            ktab(ji) = ktab(ji) + kproj(l,m)                       &
     &                     * ylm(l,m,klon(ji),klat(ji))
          ENDDO
        ENDDO

#if defined MPI
        IF (.NOT.external_vector_decomposition) THEN
          CALL MPI_ALLREDUCE (MPI_IN_PLACE, ktab, jpi,               &
     &                        MPI_DOUBLE_PRECISION,                  &
     &                        MPI_SUM,mpi_comm_sphylm,mpi_code)
        ENDIF
#endif

        END SUBROUTINE back_ylm_loc
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE regr_ylm( kregr, kwei, kobs, klon, klat, kobswei )
!----------------------------------------------------------------------
!                  ***  regr_ylm  ***
! 
! ** Purpose :   regression of observations along spherical harmonics
!
! Input: kobs : observation vector
!        klon : longitude of observations
!        klat : latitude of observations
!        kobswei : observation weight (typically: inverse obs error std)
!        kwei : weight of spherical harmonics (typically: signal std)
!
! Output: kregr : amplitude along every spherical harmonics
!
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: kobs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat, kobswei
        REAL(KIND=8), DIMENSION(0:,-jpl:), INTENT( in ) :: kwei
        REAL(KIND=8), DIMENSION(0:,-jpl:), INTENT( out ) :: kregr

        INTEGER :: jpi, k, l, m, l0, l1, jl, kmin, kmax, allocok
        INTEGER :: iter, maxiter, maxmode, overlap
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dobs
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dregr
        CHARACTER(len=80) :: rtype, texterror
        REAL(KIND=8) :: relstd, refstd, eps

        jpi=SIZE(kobs,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in regr_ylm'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in regr_ylm'
        IF (SIZE(kobswei,1).NE.jpi) STOP 'Inconsistent size in regr_ylm'
        IF (UBOUND(kregr,1).NE.jpl) STOP 'Inconsistent size in regr_ylm'
        IF (UBOUND(kregr,2).NE.jpl) STOP 'Inconsistent size in regr_ylm'
        IF (UBOUND(kwei,1).NE.jpl) STOP 'Inconsistent size in regr_ylm'
        IF (UBOUND(kwei,2).NE.jpl) STOP 'Inconsistent size in regr_ylm'

        rtype = regr_type
        maxiter = regr_maxiter
        maxmode = jpl * regr_maxbloc
        eps = regr_epsilon
        overlap = regr_overlap

        relstd = SQRT(SUM(kobs(:)*kobs(:)*kobswei(:)*kobswei(:))/jpi)
        PRINT *, 'Initial relative std',relstd

        SELECT CASE(regr_type)
        CASE('global')
! Global regression algorithm:
!   using all spherical harmonics together as basis functions

          l0 = 0
          l1 = jpl

          CALL regr_ylm_loc(kregr(l0:l1,:),kwei(l0:l1,:),kobs,klon,klat,kobswei,l0,l1)

        CASE('local')
! Local regression algorithm:
!   separate the regression into block of degrees (from l0 to l1)
!   iterate until convergence

          ALLOCATE( dobs(jpi), stat=allocok )
          IF (allocok.NE.0) STOP 'Allocation error in regr_ylm'
          ALLOCATE( dregr(0:jpl,-jpl:jpl), stat=allocok )
          IF (allocok.NE.0) STOP 'Allocation error in regr_ylm'

          kregr(:,:) = 0.0

          refstd = SQRT(SUM(kobs(:)*kobs(:)*kobswei(:)*kobswei(:))/jpi)
          DO iter=1,maxiter
            l1 = -1
            DO WHILE (l1.LT.jpl)
              l0 = MAX(0,l1+1-overlap)
              l1 = INT( SQRT(REAL((l0+1)*(l0+1)+maxmode-1,8)) - 1 )
              l1 = MIN(MAX(l1,l0+overlap),jpl)
              dregr(:,:) = 0.0
              CALL regr_ylm_loc(dregr(l0:l1,:),kwei(l0:l1,:),kobs,klon,klat,kobswei,l0,l1)
              kregr(:,:) = kregr(:,:) + dregr(:,:)
              CALL back_ylm_loc(dregr,dobs,klon,klat,l0,l1)
              kobs(:) = kobs(:) - dobs(:)
            ENDDO
            relstd = SQRT(SUM(kobs(:)*kobs(:)*kobswei(:)*kobswei(:))/jpi)
            PRINT *, 'Iteration',iter,'Relative residual std',relstd
            IF (ABS(relstd-refstd)/refstd.LE.eps) THEN
              EXIT
            ELSE
              refstd = relstd
              IF (iter.EQ.maxiter) THEN
                PRINT *, 'Warning: no convergence in regr_ylm'
              ENDIF
            ENDIF
          ENDDO

          IF (allocated(dobs)) deallocate(dobs)
          IF (allocated(dregr)) deallocate(dregr)

        CASE('sublocal')
! Sublocal regression algorithm:
!   separate the regression into block of degrees (from l0 to l1)
!   iterate until convergence before going to the next block of degrees

          ALLOCATE( dobs(jpi), stat=allocok )
          IF (allocok.NE.0) STOP 'Allocation error in regr_ylm'
          ALLOCATE( dregr(0:jpl,-jpl:jpl), stat=allocok )
          IF (allocok.NE.0) STOP 'Allocation error in regr_ylm'

          kregr(:,:) = 0.0

          jl = -1
          DO WHILE (jl.LT.jpl)
            l0 = jl + 1
            jl = INT( SQRT(REAL((l0+1)*(l0+1)+maxmode-1,8)) - 1 )
            jl = MIN(MAX(jl,l0),jl)

            refstd = SQRT(SUM(kobs(:)*kobs(:)*kobswei(:)*kobswei(:))/jpi)
            DO iter=1,maxiter
              l1 = -1
              DO WHILE (l1.LT.jl)
                l0 = MAX(0,l1+1-overlap)
                l1 = INT( SQRT(REAL((l0+1)*(l0+1)+maxmode-1,8)) - 1 )
                l1 = MIN(MAX(l1,l0+overlap),jl)
                dregr(:,:) = 0.0
                CALL regr_ylm_loc(dregr(l0:l1,:),kwei(l0:l1,:),kobs,klon,klat,kobswei,l0,l1)
                kregr(:,:) = kregr(:,:) + dregr(:,:)
                CALL back_ylm_loc(dregr,dobs,klon,klat,l0,l1)
                kobs(:) = kobs(:) - dobs(:)
              ENDDO
              relstd = SQRT(SUM(kobs(:)*kobs(:)*kobswei(:)*kobswei(:))/jpi)
              PRINT *, 'Max degree:',jl,'Iteration',iter
              PRINT *, '   Relative residual std',relstd
              IF (ABS(relstd-refstd)/refstd.LE.eps) THEN
                EXIT
              ELSE
                refstd = relstd
                IF (iter.EQ.maxiter) THEN
                  PRINT *, 'Warning: no convergence in regr_ylm'
                  PRINT *, '   at degree:',jl
                ENDIF
              ENDIF
            ENDDO

          ENDDO

          IF (allocated(dobs)) deallocate(dobs)
          IF (allocated(dregr)) deallocate(dregr)

        CASE DEFAULT
          STOP 'Bad regression type in regr_ylm'
        END SELECT

        END SUBROUTINE regr_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE regr_ylm_loc( kregr, kwei, kobs, klon, klat, kobswei, kl0, kl1 )
!----------------------------------------------------------------------
!                  ***  regr_ylm_loc  ***
! 
! ** Purpose :   regression along spherical harmonics
!                over a local range of degrees
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kobs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat, kobswei
        REAL(KIND=8), DIMENSION(kl0:,-jpl:), INTENT( in ) :: kwei
        REAL(KIND=8), DIMENSION(kl0:,-jpl:), INTENT( out ) :: kregr
        INTEGER, INTENT( in ) :: kl0, kl1

        INTEGER :: jpi, k, l, m, kmin, kmax, allocok
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vectb,vectx
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mata

        jpi=SIZE(kobs,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in regr_ylm_loc'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in regr_ylm_loc'
        IF (SIZE(kobswei,1).NE.jpi) STOP 'Inconsistent size in regr_ylm_loc'
        IF (UBOUND(kregr,1).NE.kl1) STOP 'Inconsistent size in regr_ylm_loc'
        IF (UBOUND(kregr,2).NE.jpl) STOP 'Inconsistent size in regr_ylm_loc'
        IF (UBOUND(kwei,1).NE.kl1) STOP 'Inconsistent size in regr_ylm_loc'
        IF (UBOUND(kwei,2).NE.jpl) STOP 'Inconsistent size in regr_ylm_loc'

        kmin = (kl0+1)*(kl0+1) - 1
        kmax = (kl1+1)*(kl1+1) - 1
        ALLOCATE( vectb(kmax-kmin+1), stat=allocok )
        IF (allocok.NE.0) STOP 'Allocation error in regr_ylm_loc'
        vectb(:) = 0.0_8
        ALLOCATE( vectx(kmax-kmin+1), stat=allocok )
        IF (allocok.NE.0) STOP 'Allocation error in regr_ylm_loc'
        vectx(:) = 0.0_8
        ALLOCATE( mata(kmax-kmin+1,kmax-kmin+1), stat=allocok )
        IF (allocok.NE.0) STOP 'Allocation error in regr_ylm_loc'
        mata(:,:) = 0.0_8

        CALL regr_calc_ab(mata,vectb,kwei,kobs,klon,klat,kobswei,kl0,kl1)
        CALL regr_calc_x(vectx,mata,vectb)

        DO k=kmin,kmax
          l=INT(SQRT(REAL(k,8)))
          m=k-l*l-l
          kregr(l,m) = vectx(k-kmin+1) * kwei(l,m)
        ENDDO

! --- deallocation
        IF (allocated(vectb)) deallocate(vectb)
        IF (allocated(vectx)) deallocate(vectx)
        IF (allocated(mata))  deallocate(mata)

        END SUBROUTINE regr_ylm_loc
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE regr_calc_ab( kmata, kvectb, kwei, kobs, klon, klat, kobswei, kl0, kl1 )
!----------------------------------------------------------------------
!                  ***  regr_calc_ab  ***
! 
! ** Purpose :   computation of matrix A and vector b of the linear system
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kobs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat, kobswei
        REAL(KIND=8), DIMENSION(kl0:,-jpl:), INTENT( in ) :: kwei
        REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: kmata
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: kvectb
        INTEGER, INTENT( in ) :: kl0, kl1

        INTEGER :: jpi, ji, k, k2, l, m, l2, m2, kmin, kmax
        REAL(KIND=8) :: ylm1, ylm2

        kmin = (kl0+1)*(kl0+1) - 1
        kmax = (kl1+1)*(kl1+1) - 1

        jpi=SIZE(kobs,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in regr_calc_ab'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in regr_calc_ab'
        IF (SIZE(kobswei,1).NE.jpi) STOP 'Inconsistent size in regr_calc_ab'
        IF (SIZE(kmata,1).NE.kmax-kmin+1) STOP 'Inconsistent size in regr_calc_ab'
        IF (SIZE(kmata,2).NE.kmax-kmin+1) STOP 'Inconsistent size in regr_calc_ab'
        IF (SIZE(kvectb,1).NE.kmax-kmin+1) STOP 'Inconsistent size in regr_calc_ab'
        IF (UBOUND(kwei,1).NE.kl1) STOP 'Inconsistent size in regr_calc_ab'
        IF (UBOUND(kwei,2).NE.jpl) STOP 'Inconsistent size in regr_calc_ab'

        ! Loop on spherical harmonics
        ! (reduced to one single loop for an easy parallelization)
        DO k=jproc+kmin,kmax,jpproc
          l=INT(SQRT(REAL(k,8)))
          m=k-l*l-l
          
          ! Compute A (kmata) and b (kvectb) in regression system Ax=b
          DO ji=1,jpi
            ylm1 = ylm(l,m,klon(ji),klat(ji)) * kwei(l,m) * kobswei(ji)
            kvectb(k-kmin+1) = kvectb(k-kmin+1) + kobs(ji) * kobswei(ji) * ylm1
            DO k2=kmin,k
              IF (k2.EQ.k) THEN
                ylm2 = ylm1
              ELSE
                l2=INT(SQRT(REAL(k2,8)))
                m2=k2-l2*l2-l2
                ylm2 = ylm(l2,m2,klon(ji),klat(ji)) * kwei(l2,m2) * kobswei(ji)
              ENDIF
              kmata(k-kmin+1,k2-kmin+1) = kmata(k-kmin+1,k2-kmin+1) + ylm1 * ylm2
            ENDDO
          ENDDO
        ENDDO

#if defined MPI
        CALL MPI_ALLREDUCE (MPI_IN_PLACE, kvectb, kmax,             &
     &                      MPI_DOUBLE_PRECISION,                   &
     &                      MPI_SUM,mpi_comm_sphylm,mpi_code)
        CALL MPI_ALLREDUCE (MPI_IN_PLACE, kmata, kmax*kmax,         &
     &                      MPI_DOUBLE_PRECISION,                   &
     &                      MPI_SUM,mpi_comm_sphylm,mpi_code)
#endif

!       Compute second half of symmetric matrix
        DO k=1,kmax-kmin+1
        DO k2=k+1,kmax-kmin+1
          kmata(k,k2) = kmata(k2,k)
        ENDDO
        ENDDO

!       Add identity matrix
        DO k=1,kmax-kmin+1
          kmata(k,k) = kmata(k,k) + regr_rho
        ENDDO

        END SUBROUTINE regr_calc_ab
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE regr_calc_x( kvectx, kmata, kvectb )
!----------------------------------------------------------------------
!                  ***  regr_calc_  ***
! 
! ** Purpose :   solve the linear system
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: kmata
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kvectb
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: kvectx

        INTEGER :: jk, jpk
        REAL(KIND=8) :: sc

        jpk=SIZE(kvectx,1)
        IF (SIZE(kvectb,1).NE.jpk) STOP 'Inconsistent size in regr_calc_x'
        IF (SIZE(kmata,1).NE.jpk) STOP 'Inconsistent size in regr_calc_x'
        IF (SIZE(kmata,2).NE.jpk) STOP 'Inconsistent size in regr_calc_x'

        sc=maxval(abs(kmata))
        kmata=kmata/sc
        CALL invmat(kmata)
        DO jk=1,jpk
          kvectx(jk)=DOT_PRODUCT(kmata(jk,:),kvectb(:))
        ENDDO
        kvectx=kvectx/sc

        END SUBROUTINE regr_calc_x
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE disp_ylm( ktab, klon, klat, kl, km )
!----------------------------------------------------------------------
!                  ***  init_ylm  ***
! 
! ** Purpose :   output one single spherical harmonics
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: klon, klat
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ktab
        INTEGER, INTENT( in ) :: kl, km

        INTEGER :: jpi, ji, k

        jpi=SIZE(ktab,1)
        IF (SIZE(klon,1).NE.jpi) STOP 'Inconsistent size in disp_ylm'
        IF (SIZE(klat,1).NE.jpi) STOP 'Inconsistent size in disp_ylm'
        IF (ABS(km).GT.kl) STOP 'Inconsistent order in disp_ylm'

        ktab(:)=0.0
        DO ji=1,jpi
          ktab(ji) =  ylm(kl,km,klon(ji),klat(ji))
        ENDDO

        END SUBROUTINE disp_ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION ylm(kl,km,klon,klat)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: kl, km
        REAL(KIND=8), INTENT(in) :: klon, klat

        REAL(KIND=8) :: ylm, xlon, xlat
        INTEGER :: jlat

        xlon = klon * deg2rad
        xlat = (klat-lat0)/dlat

        jlat = INT(xlat)
        xlat = xlat - jlat

        ylm = pleg(jlat,kl,ABS(km))
        IF (jlat < jplat) ylm = ylm + xlat * ( pleg(jlat+1,kl,ABS(km)) &
     &                                       - pleg(jlat,kl,ABS(km)) )
        IF (km>0) ylm = ylm * COS( km * xlon ) * SQRT(2.0)
        IF (km<0) ylm = ylm * SIN( km * xlon ) * SQRT(2.0)

        END FUNCTION ylm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE plm(x,k)
!----------------------------------------------------------------------
!                  ***  SUBROUTINE  plm ***
! 
! ** Purpose :   Compute associated Legendre polynomails
!                of degrees up to jpl, for argument x
!                (assuming that |x|<=1)
!
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT(in) :: x
        INTEGER, INTENT( in ) :: k

        INTEGER :: l ! degree of polynomial
        INTEGER :: m ! order of polynomial
        INTEGER :: i
        REAL(KIND=8) :: xq, fac, a, b

        ! Set all polynomials to zero
        pleg(k,:,:) = 0.0

        ! Degree l=0 and order m=0
        pleg(k,0,0) = 1.0

        ! x = cosine, xq = sine of theta
        xq = SQRT( 1.0 - x*x )

        ! Recursive computation for m=l (order equal to degree)
        DO l = 1, jpl
          pleg(k,l,l) = SQRT(REAL(2*l+1,8)/ REAL(2*l,8)) &
     &                  * xq * pleg(k,l-1,l-1)
        ENDDO

        ! Recursive computation for m<l (other orders)
        DO m = 0, jpl
        DO l = m+1, jpl
          a = REAL(2*l-1,8) * REAL(2*l+1,8)
          a = a / (REAL(l-m,8) * REAL(l+m,8))
          pleg(k,l,m) = x * SQRT(a) * pleg(k,l-1,m)
          IF (l.GT.m+1) THEN
            b = REAL(2*l+1,8) * REAL(l+m-1,8) * REAL(l-m-1,8)
            b = b / (REAL(l-m,8) * REAL(l+m,8) * REAL(2*l-3,8))
            pleg(k,l,m) = pleg(k,l,m) - SQRT(b) * pleg(k,l-2,m)
          ENDIF
        ENDDO
        ENDDO

        END SUBROUTINE plm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE invmat(mat)
!----------------------------------------------------------------------
!                  ***  SUBROUTINE invmat ***
!
! Compute inverse of a symmetric matrix C using Choleski decomposition
!
! Argument:  mat : C (input and output)
!
!----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=8), dimension(:,:), intent(inout) :: mat
      REAL(KIND=8), dimension(:,:), allocatable :: invl
      REAL(KIND=8), dimension(:), allocatable :: p
      REAL(KIND=8) :: sum
      INTEGER :: jprsize,jr1,jr2,jk
      INTEGER :: allocok

      jprsize = size(mat,1)
      IF (jprsize.NE.size(mat,2)) STOP 'Inconsistent size in invmat'

      allocate(invl(jprsize,jprsize),stat=allocok)
      IF (allocok.NE.0) STOP 'Allocation error in invmat'
      allocate(p(jprsize),stat=allocok)
      IF (allocok.NE.0) STOP 'Allocation error in invmat'

! Choleski decomposision
! ----------------------
! Compute L such that: C = L L^T
      DO jr1=1,jprsize
      DO jr2=1,jprsize
        sum = mat(jr1,jr2)
        DO jk=jr1-1,1,-1
          sum = sum - mat(jr1,jk) * mat(jr2,jk)
        ENDDO
        IF (jr1.EQ.jr2) THEN
          IF (sum.LE.0.) STOP 'Error in Choleski decomposition'
          p(jr1) = SQRT(sum)
        ELSE
          mat(jr2,jr1)=sum/p(jr1)
        ENDIF
      ENDDO
      ENDDO

      DO jr2=1,jprsize
         mat(1:jr2-1,jr2) = 0.0
      ENDDO

! Inverse computation
! -------------------
      ! Inverse L
      DO jr1=1,jprsize
        mat(jr1,jr1)=1.0/p(jr1)
        DO jr2=jr1+1,jprsize
          sum = 0.0
          DO jk=jr1,jr2-1
            sum = sum - mat(jr2,jk) * mat(jk,jr1)
          ENDDO
          mat(jr2,jr1) = sum / p(jr2)
        ENDDO
      ENDDO
      invl = mat

      ! Inverse A
      DO jr1=1,jprsize
      DO jr2=1,jprsize
        mat(jr1,jr2) = DOT_PRODUCT(invl(:,jr1),invl(:,jr2))
      ENDDO
      ENDDO

      IF (allocated(invl)) deallocate(invl)
      IF (allocated(p)) deallocate(p)

      END SUBROUTINE invmat
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_sphylm
