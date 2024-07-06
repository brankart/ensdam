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
!                        MODULE SCORE_CRPS
!
!---------------------------------------------------------------------
! Computation of CRPS score by comparing ensemble simulation to verification data
! by Jean-Michel Brankart, November 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! crps_score : compute CRPS score (with option to partition the data)
! crps_cumul : accumulate data to prepare the final computation of the score
! crps_final : compute final score from accumulated data
! ----------------------------------------------------------------------
MODULE ensdam_score_crps
#ifdef MPI_MODULE
      use mpi
#endif
      IMPLICIT NONE
      PRIVATE
#ifdef MPI_INCLUDE
      include "mpif.h"
#endif

      PUBLIC crps_score, crps_cumul, crps_final

      INTERFACE crps_score
        MODULE PROCEDURE crps_score_global, crps_score_partition
      END INTERFACE

      REAL(KIND=8), PUBLIC, SAVE  :: crps_missing_value = -9999.

      ! Public definitions needed by python/julia APIs
      PUBLIC crps_score_global, crps_score_partition

      ! Definition for MPI
#if defined MPI
      INTEGER, PUBLIC, SAVE  :: mpi_comm_score_crps=mpi_comm_world   ! definition of module global communicator
      INTEGER, SAVE :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE crps_score_global( crps, reliability, resolution, ens, verif )
!----------------------------------------------------------------------
! ** Purpose :   compute CRPS score
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to verification data)
!         verif  : verification data
!         crps        : CRPS score for each region (score = reliability + resolution)
!         reliability : reliability part of CRPS
!         resolution  : resolution part of CRPS
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( out ) :: crps
      REAL(KIND=8), INTENT( out ) :: reliability
      REAL(KIND=8), INTENT( out ) :: resolution
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif

      INTEGER :: jpi,jpm,nbr,ji,allocstat
      REAL(KIND=8), DIMENSION(:), allocatable :: aa,bb

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in score_crps'

      allocate( aa(0:jpm), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'
      allocate( bb(0:jpm), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'

      aa = 0. ; bb = 0. ; nbr = 0

      DO ji=1,jpi
        nbr = nbr + 1
        CALL crps_cumul(ens(ji,:),verif(ji),aa(0:),bb(0:))
      ENDDO

#if defined MPI
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, aa, jpm+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, bb, jpm+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, 1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
#endif

      IF (nbr.GT.0) THEN
        aa(0:) = aa(0:) / nbr
        bb(0:) = bb(0:) / nbr
        CALL crps_final(aa(0:),bb(0:),reliability,resolution,crps)
      ENDIF

      deallocate(aa,bb)

      END SUBROUTINE crps_score_global
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE crps_score_partition( crps, reliability, resolution, ens, verif, partition )
!----------------------------------------------------------------------
! ** Purpose :   compute CRPS score with partition of the data
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to verification data)
!         verif  : verification data
!         crps        : CRPS score for each region (score = reliability + resolution)
!         reliability : reliability part of CRPS
!         resolution  : resolution part of CRPS
!         partition   : partition of verification data
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: crps
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: reliability
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: resolution
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif
      INTEGER, DIMENSION(:), INTENT( in ) :: partition

      INTEGER :: jpi,jpm,submin,submax,jsub,ji,allocstat
      REAL(KIND=8), DIMENSION(:,:), allocatable :: aa,bb
      INTEGER, DIMENSION(:), allocatable :: nbr

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in score_crps'
      IF (SIZE(partition,1).NE.jpi) STOP 'Inconsistent size in score_crps'

      submin = MINVAL(partition)  ! Maximum index of subdomains in verification data
      submax = MAXVAL(partition)  ! Maximum index of subdomains in verification data

      IF (LBOUND(crps,1).NE.submin) STOP 'Inconsistent array bounds in score_crps'
      IF (UBOUND(crps,1).NE.submax) STOP 'Inconsistent array bounds in score_crps'
      IF (LBOUND(reliability,1).NE.submin) STOP 'Inconsistent array bounds in score_crps'
      IF (UBOUND(reliability,1).NE.submax) STOP 'Inconsistent array bounds in score_crps'
      IF (LBOUND(resolution,1).NE.submin) STOP 'Inconsistent array bounds in score_crps'
      IF (UBOUND(resolution,1).NE.submax) STOP 'Inconsistent array bounds in score_crps'

      allocate( aa(0:jpm,submin:submax), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'
      allocate( bb(0:jpm,submin:submax), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'
      allocate( nbr(submin:submax), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'

      aa = 0. ; bb = 0. ; nbr = 0

      DO ji=1,jpi
        nbr(partition(ji)) = nbr(partition(ji)) + 1
        CALL crps_cumul(ens(ji,:),verif(ji),aa(0:,partition(ji)),bb(0:,partition(ji)))
      ENDDO

#if defined MPI
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, aa, (jpm+1)*(submax-submin+1), MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, bb, (jpm+1)*(submax-submin+1), MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, submax-submin+1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_crps,mpi_code)
#endif

      DO jsub=submin,submax
        IF (nbr(jsub).GT.0) THEN
          aa(0:,jsub) = aa(0:,jsub) / nbr(jsub)
          bb(0:,jsub) = bb(0:,jsub) / nbr(jsub)
          CALL crps_final(aa(0:,jsub),bb(0:,jsub),reliability(jsub),resolution(jsub),crps(jsub))
        ELSE
          reliability(jsub) = crps_missing_value
          resolution(jsub) = crps_missing_value
          crps(jsub) = crps_missing_value
        ENDIF
      ENDDO

      deallocate(aa,bb,nbr)

      END SUBROUTINE crps_score_partition
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE crps_cumul(ens,a,aa,bb)
!----------------------------------------------------------------------
! ** Purpose :   Accumulate one more piece of information to compute CRPS score
! 
! ** Arguments :
!         e  : input ensemble
!         a  : verification
!         aa : array to update
!         bb : array to update
!----------------------------------------------------------------------
      use ensdam_anaqua
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ens
      REAL(KIND=8), INTENT(IN) :: a
      REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: aa, bb

      INTEGER :: n, i, allocstat
      REAL(KIND=8), DIMENSION(:), allocatable :: e

      n = size(ens)

      allocate( e(n), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_crps'
      e = ens

! Sort input ensemble
      CALL heapsort(e)

! Verfication smaller than all ensemble members
      IF(a.LT.e(1)) bb(0)=bb(0)+1.0
      IF(a.LT.e(1)) aa(0)=aa(0)+(e(1)-a)

! Verification inside ensemble range
      DO i=1,n-1
        bb(i)=bb(i)+MAX(MIN(a,e(i+1))-e(i),0._8)
        aa(i)=aa(i)+MAX(e(i+1)-MAX(a,e(i)),0._8)
      ENDDO

! Verfication larger than all ensemble members
      IF(a.GT.e(n)) bb(n)=bb(n)+(a-e(n))
      IF(a.GT.e(n)) aa(n)=aa(n)+1.0

      deallocate(e)

      END SUBROUTINE crps_cumul
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE crps_final(aa,bb,reli,resol,crps)
!----------------------------------------------------------------------
! ** Purpose : Compute CRPS score from previously accumulated information
! 
! ** Arguments :
!         aa : accumulated information
!         bb : accumulated information
!         reli  : reliability score
!         resol : resolution score
!         crps  : CRPS score (= reli + resol)
!----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: aa, bb
      REAL(KIND=8), INTENT(OUT) :: reli, resol, crps

      INTEGER :: n, i
      REAL(KIND=8) :: gi, oi, p

      n = size(aa) - 1

      crps = 0. ; reli = 0. ; resol = 0.

! Reinterpretation of bb(0) and aa(n)
! (note that this does not contribute to total CRPS)
      IF (bb(0).NE.0.) bb(0)=aa(0)*(1./bb(0)-1.)
      IF (aa(n).NE.0.) aa(n)=bb(n)*(1./aa(n)-1.)

! Compute components of CRPS
      DO i=0,n
        gi = bb(i) + aa(i)
        IF (gi.NE.0.) oi = aa(i)/gi
        p = REAL(i,8)/REAL(n,8)

        crps  = crps  + bb(i)*p*p+aa(i)*(1.-p)*(1.-p)
        reli  = reli  + gi*(oi-p)**2
        resol = resol + gi*oi*(1.-oi)
      ENDDO

      END SUBROUTINE crps_final
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_score_crps
