! ----------------------------------------------------------------------
! Copyright: CNRS - Universite de Grenoble Alpes
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
!                        MODULE STOTGE
!
!---------------------------------------------------------------------
! Module to sample truncated Gaussian or truncated exponential distributions
! by Jean-Michel Brankart, June 2011,
! revised September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ranv_tg  : Multidimensional truncated Gaussian random vector generator
! ran_tg   : One-dimensional truncated Gaussian random number generator
! ran_te   : One-dimensional truncated exponential random number generator
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE ensdam_stotge
   use ensdam_storng

   IMPLICIT NONE
   PRIVATE

! Public functions/subroutines
   PUBLIC :: ranv_tg, ran_tg, ran_te

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ranv_tg(tgvsmpl,matArm,vecbm)
!CC---------------------------------------------------------------------
!CC
!CC  Purpose : Multidimensional truncated Gaussian random vector generator
!CC  -------   (sample p(x)  ~  exp (- xT x / 2 ) , A x <= b)
!CC
!CC  Method :  Gibbs sampler, using one-dimensional truncated Gaussian
!CC  ------    random number generator
!CC
!CC  Input :   vecbm   : b vector (jpmsize)           ; jpmsize = nbr of constraints
!CC  -----     matArm  : A matrix (jprsize x jpmsize) ; jprsize = size of x vector
!CC
!CC  Output :  tgvsmpl : x sample (jpsmpl  x jprsize) ; jpsmpl  = size of the sample
!CC  ------
!CC---------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL(KIND=8), dimension(:,:), intent(inout) :: tgvsmpl
      REAL(KIND=8), dimension(:,:), intent(in) :: matArm
      REAL(KIND=8), dimension(:), intent(in) :: vecbm
!
      REAL(KIND=8), dimension(1) :: tgran
      REAL(KIND=8) :: a, b, abnew
      INTEGER :: jprsize, jpmsize, jpsmpl, jsmpl, jr, jm, jr1
!
      jpsmpl=size(tgvsmpl,1)
      jprsize=size(tgvsmpl,2)
      jpmsize=size(matArm,2)
      IF (jprsize.NE.size(matArm,1)) STOP 'ranv_tg: Inconsistent vector size'
      IF (jpmsize.NE.size(vecbm,1)) STOP 'ranv_tg: Inconsistent vector size'
!
      DO jsmpl=2,jpsmpl
        jr=MOD(jsmpl-2,jprsize)+1
!
! Compute interval bound for the jr conditional tg pdf
! by iterating on the inequality constraints
        a=-huge(1.0) ; b=huge(1.0)
        DO jm=1,jpmsize
! 
          abnew=DOT_PRODUCT(matArm(:,jm),tgvsmpl(jsmpl-1,:))
          abnew=abnew-matArm(jr,jm)*tgvsmpl(jsmpl-1,jr)
          abnew=vecbm(jm)-abnew
!
          IF (matArm(jr,jm).LT.0.0) THEN
            abnew=abnew/matArm(jr,jm)
            IF (abnew.GT.a) a=abnew
          ELSEIF (matArm(jr,jm).GT.0.0) THEN
            abnew=abnew/matArm(jr,jm)
            IF (abnew.LT.b) b=abnew
          ELSE
            IF (abnew.LT.0.0) STOP 'ranv_tg: Impossible initial vector'
          ENDIF
!
        ENDDO
        IF (b.LT.a) STOP 'ranv_tg: Impossible inequality constraints'
!
! Sample 1D conditional tg pdf
        CALL ran_tg(tgran,a,b)
!
! Build next draw of the Gibbs sample
        tgvsmpl(jsmpl,1:jprsize)=tgvsmpl(jsmpl-1,1:jprsize)
        tgvsmpl(jsmpl,jr)=tgran(1)
      ENDDO
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ran_tg(tgsmpl,aa,bb)
!
! One-dimensional truncated Gaussian random number generator
! (sample p(y)  ~  exp(-y*y) , aa < y < bb)
! The algorithm is taken from Geweke, 1991.
!
      IMPLICIT NONE
! --- Variable declaration 
      REAL(KIND=8), dimension(:), intent(inout) :: tgsmpl
      REAL(KIND=8), intent(in) :: aa, bb
!
      INTEGER :: jpsmpl, jsmpl
      REAL(KIND=8), parameter :: pi=3.1415926535897932384526
      REAL(KIND=8), parameter :: t1=0.15, t2=2.18, t3=0.725, t4=0.45
      REAL(KIND=8) :: a, b, g0, teran, ratio, uranab, guranab
      REAL(KIND=8) :: uran, gran
      LOGICAL :: switch, found
!
      g0=SQRT(2.0*pi)
      jpsmpl=size(tgsmpl,1)
!
! Check the validity of interval bounds
      IF (aa.GT.bb) STOP 'ran_tg: empty interval'
      IF (aa.EQ.bb) THEN
        tgsmpl(1:jpsmpl)=aa ; RETURN
      ENDIF
!
! Use symmetry to ensure |a|<|b| and a<b
      IF (ABS(aa).GT.ABS(bb)) THEN
        a=-bb ; b=-aa ; switch=.TRUE.
      ELSE
        a=aa  ; b=bb  ; switch=.FALSE.
      ENDIF
!
      IF (b.EQ.huge(1.0)) THEN
! sample in the interval [a,infinity]
! -----------------------------------
        IF (a.EQ.-huge(1.0)) THEN
! normal sampling (if a=-infinity)
!         print *, 'normal sampling (a=-inf, b=inf)'
          DO jsmpl=1,jpsmpl
            CALL kiss_gaussian(gran)
            tgsmpl(jsmpl)=gran
          ENDDO
        ELSEIF (a.LE.t4) THEN
          IF (a.LT.0.0) THEN
! normal rejection sampling (if a<0)
!           print *, 'normal rejection sampling (a<0, b=inf)'
            DO jsmpl=1,jpsmpl
              found=.FALSE.
              DO WHILE (.NOT.found)
                CALL kiss_gaussian(gran)
                IF (gran.GE.a) THEN
                  tgsmpl(jsmpl)=gran ; found=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ELSE
! half normal rejection sampling (if 0<a<t4)
!           print *, 'half normal rejection sampling (0<a<t4, b=inf)'
            DO jsmpl=1,jpsmpl
              found=.FALSE.
              DO WHILE (.NOT.found)
                CALL kiss_gaussian(gran)
                IF (abs(gran).GE.a) THEN
                  tgsmpl(jsmpl)=abs(gran) ; found=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ELSE
! exponential rejection sampling (if a>t4)
!         print *, 'exponential rejection sampling (a>t4, b=inf)'
          DO jsmpl=1,jpsmpl
            found=.FALSE.
            DO WHILE (.NOT.found)
              CALL ran_te(teran,a)
              CALL kiss_uniform(uran)
              ratio=EXP(-(teran-a)*(teran-a)/2.0)
              IF (uran.LE.ratio) THEN
                tgsmpl(jsmpl)=teran ; found=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ELSEIF (a*b.LT.0.0) THEN
! sample in the interval [a,b] if a<0<b
! -------------------------------------
        IF ((-a*a.LE.LOG(t1*g0)).OR.(-b*b.LE.LOG(t1*g0))) THEN
! normal rejection sampling (if g(a)<t1 or g(b)<t1)
!         print *, 'normal rejection sampling (a<0<b, ga<t1 or gb<t1)'
          DO jsmpl=1,jpsmpl
            found=.FALSE.
            DO WHILE (.NOT.found)
              CALL kiss_gaussian(gran)
              IF ((gran.GE.a).AND.(gran.LE.b)) THEN
                tgsmpl(jsmpl)=gran ; found=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ELSE
! uniform rejection sampling (if g(a)>t1 and g(b)>t1)
!         print *, 'uniform rejection sampling (a<0<b, ga>t1 and gb>t1)'
          DO jsmpl=1,jpsmpl
            found=.FALSE.
            DO WHILE (.NOT.found)
              CALL kiss_uniform(uran)
              uranab=a+(b-a)*uran
              guranab=EXP(-uranab*uranab)
              CALL kiss_uniform(uran)
              IF (uran.LE.guranab) THEN
                tgsmpl(jsmpl)=uranab ; found=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ELSE
! sample in the interval [a,b] if 0<a<b
! -------------------------------------
        IF (b*b-a*a.LE.LOG(t2)) THEN
! uniform rejection sampling (if g(a)/g(b) < t2)
!         print *, 'uniform rejection sampling (0<a<b, ga/gb<t2)'
          DO jsmpl=1,jpsmpl
            found=.FALSE.
            DO WHILE (.NOT.found)
              CALL kiss_uniform(uran)
              uranab=a+(b-a)*uran
              guranab=EXP(a*a-uranab*uranab)
              CALL kiss_uniform(uran)
              IF (uran.LE.guranab) THEN
                tgsmpl(jsmpl)=uranab ; found=.TRUE.
              ENDIF
            ENDDO
          ENDDO
        ELSE
          IF (a.LE.t3) THEN
! half normal rejection sampling
!           print *, 'half normal rejection sampling (0<a<t3, ga/gb>t2)'
            DO jsmpl=1,jpsmpl
              found=.FALSE.
              DO WHILE (.NOT.found)
                CALL kiss_gaussian(gran)
                IF ((abs(gran).GE.a).AND.(abs(gran).LE.b)) THEN
                  tgsmpl(jsmpl)=abs(gran) ; found=.TRUE.
                ENDIF
              ENDDO
            ENDDO
          ELSE
! exponential rejection sampling
!           print *, 'exponential rejection sampling (t3<a<b, ga/gb>t2)'
            DO jsmpl=1,jpsmpl
              found=.FALSE.
              DO WHILE (.NOT.found)
                CALL ran_te(teran,a)
                CALL kiss_uniform(uran)
                ratio=EXP(-(teran-a)*(teran-a)/2.0)
                IF (uran.LE.ratio) THEN
                  IF (teran.LE.b) THEN
                    tgsmpl(jsmpl)=teran ; found=.TRUE.
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!
      IF (switch) tgsmpl(:)=-tgsmpl(:)

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ran_te(teran,a)
!
! Truncated exponential random number generator
! (sample p(y) = a * exp(a*a) * exp(-ay) , y > a)
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration 
      REAL(KIND=8), intent(out) :: teran
      REAL(KIND=8), intent(in) :: a

      REAL(KIND=8) :: uran

      CALL kiss_uniform(uran)

      teran=a-LOG(uran)/a

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE ensdam_stotge
