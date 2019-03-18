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
!                        MODULE ANATRA
!
!---------------------------------------------------------------------
! Perform anamorphosis transformation:
! -Forward and Backward transformation
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ana_forward : forward anamorphosis transformation
! ana_backward : backward anamorphosis transformation
! ----------------------------------------------------------------------
MODULE ensdam_anatra
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC ana_forward, ana_backward

      INTERFACE ana_forward
        MODULE PROCEDURE ana_forward_ensemble, ana_forward_vector, ana_forward_variable
      END INTERFACE

      INTERFACE ana_backward
        MODULE PROCEDURE ana_backward_ensemble, ana_backward_vector, ana_backward_variable
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC ana_forward_ensemble, ana_forward_vector, ana_forward_variable
      PUBLIC ana_backward_ensemble, ana_backward_vector, ana_backward_variable

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_forward_ensemble( ens, qua, quaref, rank )
!----------------------------------------------------------------------
!                  ***  ana_forward_ensemble  ***
! 
! ** Purpose :   forward anamorphosis transformation of input ensemble
! 
! ** Arguments :
!         ens      : input ensemble to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!         rank     : rank to define each ensemble member inside probability concentration
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: ens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: rank

        INTEGER :: jpi,jpm,jpq,ji,jm,jq
        REAL(KIND=8) :: rr

        jpi = SIZE(ens,1)
        jpm = SIZE(ens,2)
        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpi) STOP 'Inconsistent size in ana_forward'
        IF (SIZE(qua,2).NE.jpq) STOP 'Inconsistent size in ana_forward'

        DO jm = 1,jpm
          IF (present(rank)) THEN
            CALL ana_forward_vector(ens(:,jm),qua,quaref,rank=rank(jm))
          ELSE
            rr = ( jm - 0.5 ) / jpm
            CALL ana_forward_vector(ens(:,jm),qua,quaref,rank=rr)
          ENDIF
        ENDDO

        END SUBROUTINE ana_forward_ensemble
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_forward_vector( vct, qua, quaref, rank )
!----------------------------------------------------------------------
!                  ***  ana_forward_vector  ***
! 
! ** Purpose :   forward anamorphosis transformation of input vector
! 
! ** Arguments :
!         vct      : input vector to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!         rank     : rank to define transformed vector inside probability concentration
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: vct
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref
        REAL(KIND=8), INTENT( in ), OPTIONAL :: rank

        INTEGER :: jpi,jpq,ji,jq

        jpi = SIZE(vct,1)
        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpi) STOP 'Inconsistent size in ana_forward'
        IF (SIZE(qua,2).NE.jpq) STOP 'Inconsistent size in ana_forward'

        DO ji = 1,jpi
          IF (present(rank)) THEN
            CALL ana_forward_variable(vct(ji),qua(ji,:),quaref,rank=rank)
          ELSE
            CALL ana_forward_variable(vct(ji),qua(ji,:),quaref)
          ENDIF
        ENDDO

        END SUBROUTINE ana_forward_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_forward_variable( var, qua, quaref, rank )
!----------------------------------------------------------------------
!                  ***  ana_forward_variable  ***
! 
! ** Purpose :   forward anamorphosis transformation of input variable
! 
! ** Arguments :
!         var      : input variable to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!         rank     : rank to define transformed variable inside probability concentration
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( inout ) :: var
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref
        REAL(KIND=8), INTENT( in ), OPTIONAL :: rank

        INTEGER :: jpq,jq,jqinf,jqsup
        REAL(KIND=8) :: a, b, rq, rr

        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpq) STOP 'Inconsistent size in ana_forward'

        IF (present(rank)) THEN
          ! use rank provided in argument
          rr = rank
        ELSE
          ! default rank is random
          call kiss_uniform(rr)
        ENDIF

        ! Deal with possible special values
        IF (var.eq.huge(var)) THEN
          return
        ELSE
          call qualoc_fwd(jqinf,jqsup,var,qua(:))
        ENDIF

        IF (jqinf.EQ.jqsup) THEN
          ! Non-empty interval between quantiles

          jq = jqsup
          IF (jq.LE.0) THEN
            ! Below lowest quantile
            var=quaref(1)
          ELSEIF (jq.GE.jpq) THEN
            ! Above highest quantile
            var=quaref(jpq)
          ELSE
            a = (var-qua(jq))/(qua(jq+1)-qua(jq))
            var=quaref(jq)+a*(quaref(jq+1)-quaref(jq))
          ENDIF

        ELSE
          ! Empty interval between ensemble quantiles
          ! This deals with discontinuous cdf
          ! ens(ji,jm) = qua(ji,jqinf) = qua(ji,jqsup)

          IF (jqinf.LE.0) jqinf=1
          IF (jqsup.LE.0) jqsup=1
          IF (jqinf.GT.jpq) jqinf=jpq
          IF (jqsup.GT.jpq) jqsup=jpq

          ! Transform to the specified rank in the interval of quantiles
          rq = real(jqinf) + rr*real(jqsup-jqinf)
          jq = int(rq)

          b = rq - jq
          a  = 1.0 - b

          var=a*quaref(jq)
          IF (b.GT.0.) THEN
            var=var+b*quaref(jq+1)
          ENDIF

        ENDIF

        END SUBROUTINE ana_forward_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_backward_ensemble( ens, qua, quaref )
!----------------------------------------------------------------------
!                  ***  ana_forward  ***
! 
! ** Purpose :   backward anamorphosis transformation of input ensemble
! 
! ** Arguments :
!         ens      : ensemble to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: ens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref

        INTEGER :: jpi,jpm,jpq,ji,jm,jq

        jpi = SIZE(ens,1)
        jpm = SIZE(ens,2)
        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpi) STOP 'Inconsistent size in ana_backward'
        IF (SIZE(qua,2).NE.jpq) STOP 'Inconsistent size in ana_backward'

        DO jm = 1,jpm
          CALL ana_backward_vector( ens(:,jm), qua, quaref )
        ENDDO

        END SUBROUTINE ana_backward_ensemble
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_backward_vector( vct, qua, quaref )
!----------------------------------------------------------------------
!                  ***  ana_forward_vector  ***
! 
! ** Purpose :   backward anamorphosis transformation of input vector
! 
! ** Arguments :
!         vct      : vector to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: vct
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref

        INTEGER :: jpi,jpq,ji,jq

        jpi = SIZE(vct,1)
        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpi) STOP 'Inconsistent size in ana_backward'
        IF (SIZE(qua,2).NE.jpq) STOP 'Inconsistent size in ana_backward'

        DO ji = 1,jpi
          CALL ana_backward_variable( vct(ji), qua(ji,:), quaref )
        ENDDO

        END SUBROUTINE ana_backward_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_backward_variable( var, qua, quaref )
!----------------------------------------------------------------------
!                  ***  ana_forward_variable  ***
! 
! ** Purpose :   backward anamorphosis transformation of input variable
! 
! ** Arguments :
!         var      : variable to be transformed
!         qua      : ensemble quantiles (defining the piecewise linear transformation)
!         quaref   : corresponding quantiles of the reference/target distribution
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( inout ) :: var
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref

        INTEGER :: jpq,jq,jqinf,jqsup
        REAL(KIND=8) :: a, rq

        jpq = SIZE(quaref,1)
        IF (SIZE(qua,1).NE.jpq) STOP 'Inconsistent size in ana_backward'

        call qualoc_backwd(jq,var,quaref(:))

        IF (jq.LE.0) THEN
          ! Below lowest quantile
          jq = 1 ; a = 0.
        ELSEIF (jq.GE.jpq) THEN
          ! Above highest quantile
          jq = jpq-1 ; a = 1.
        ELSE
          ! Quantiles of the reference cdf are assumed always distinct
          a = (var-quaref(jq))/(quaref(jq+1)-quaref(jq))
        ENDIF

        var=qua(jq)+a*(qua(jq+1)-qua(jq))

        END SUBROUTINE ana_backward_variable
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE qualoc_fwd(jpinf,jpsup,v,qua)
!CC---------------------------------------------------------------------
!CC
!CC  Purpose : Localize value in the ordered list of quantiles;
!CC            if the quantiles are distinct: jpsup = jpinf
!CC                      = index of largest quantiles smallest or equal to v
!CC            if not, jpsup = index of largest quantile equal to v,
!CC                    jpinf = idex of smallest quantile equal to v,
!CC
!CC  Method : bisection in the sequence of quantile
!CC  ------
!CC
!CC  Input : v    : value to localize
!CC  -----   qua  : sequence of quantiles
!CC
!CC---------------------------------------------------------------------
      IMPLICIT NONE

      REAL(kind=8), intent(in) :: v
      REAL(kind=8), dimension(:), intent(in) :: qua
      INTEGER, intent(out) :: jpinf,jpsup

      INTEGER :: jp,jp0,jp1,jperc,jpp

      jpp=size(qua,1)

      jp0=0 ; jp1=jpp+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GE.qua(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO
      jpsup=jp

      jp0=0 ; jp1=jpsup ; jp=(jp0+jp1+1)/2
      DO WHILE (jp1.NE.jp)
        IF (v.GT.qua(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1+1)/2
      ENDDO
      jpinf=jp

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE qualoc_backwd(jp,v,qua)
!CC---------------------------------------------------------------------
!CC
!CC  Purpose : Localize value in the sequence of quantile
!CC            (quantiles are here assumed distinct)
!CC
!CC  Method : bisection in the sequence of quantile
!CC  ------
!CC
!CC  Input : v    : value to localize
!CC  -----   qua  : sequence of quantiles
!CC
!CC---------------------------------------------------------------------
      IMPLICIT NONE

      REAL(kind=8), intent(in) :: v
      REAL(kind=8), dimension(:), intent(in) :: qua
      INTEGER, intent(out) :: jp

      INTEGER :: jp0,jp1,jpp

      jpp=size(qua,1)

      jp0=0 ; jp1=jpp+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GE.qua(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE ensdam_anatra
