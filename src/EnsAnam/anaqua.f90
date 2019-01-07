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
!                        MODULE ANAQUA
!
!---------------------------------------------------------------------
! Computation of ensemble quantiles
! to perform anamorphosis transformation
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ens_quantiles : compute ensemble quantiles
! ----------------------------------------------------------------------
MODULE ensdam_anaqua
      IMPLICIT NONE
      PRIVATE

      PUBLIC ens_quantiles, heapsort, heapsort2

      INTERFACE ens_quantiles
        MODULE PROCEDURE ens_quantiles_vector, ens_quantiles_variable
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC ens_quantiles_vector, ens_quantiles_variable

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ens_quantiles_vector( qua, ens, quadef, enswei, ensweiloc )
!----------------------------------------------------------------------
!                  ***  ens_quantiles  ***
! 
! ** Purpose :   compute ensemble quantiles
!
! ** Arguments :
!         qua      : ensemble quantiles
!         ens      : input ensemble
!         quadef   : definition of the quantiles used for anamorphosis
!         enswei   : weight of ensemble members (optional)
!                    (default: all members have the same weight)
!         ensweiloc : local weight of ensemble members (optional)
!                     (default: all members have the same weight everywhere)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: qua
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quadef
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: enswei
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ), OPTIONAL :: ensweiloc

        REAL(KIND=8), DIMENSION(:), allocatable :: ens_sort, enswei_sort, enswei_ave
        INTEGER :: jpi,jpm,jpq,ji,jm,jq,allocstat

        jpi = SIZE(ens,1)
        jpm = SIZE(ens,2)
        jpq = SIZE(quadef,1)
        IF (SIZE(qua,1).NE.jpi) STOP 'Inconsistent size in ens_quantiles'
        IF (SIZE(qua,2).NE.jpq) STOP 'Inconsistent size in ens_quantiles'
        IF (present(enswei)) THEN
          IF (SIZE(enswei,1).NE.jpm) STOP 'Inconsistent size in ens_quantiles'
        ENDIF
        IF (present(ensweiloc)) THEN
          IF (SIZE(ensweiloc,1).NE.jpi) STOP 'Inconsistent size in ens_quantiles'
          IF (SIZE(ensweiloc,2).NE.jpm) STOP 'Inconsistent size in ens_quantiles'
        ENDIF

        IF (MINVAL(quadef).LT.0.0) STOP 'Bad quantile definition in ens_quantiles'
        IF (MAXVAL(quadef).GT.1.0) STOP 'Bad quantile definition in ens_quantiles'

        allocate ( ens_sort(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'

        IF ((.NOT.present(enswei)).AND.(.NOT.present(ensweiloc))) THEN

          DO ji = 1,jpi

            ! Sort input ensemble
            ens_sort(:) = ens(ji,:)
            CALL heapsort(ens_sort)

            ! Get quantiles
            CALL get_quantiles( qua(ji,:), ens_sort, quadef )

          ENDDO

        ELSE

          allocate ( enswei_sort(1:jpm), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'
          allocate ( enswei_ave(1:jpm), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'

          DO ji = 1,jpi

            ! Get appropriate weights
            IF (present(ensweiloc)) THEN
              enswei_sort(:) = ensweiloc(ji,:)
            ELSE
              enswei_sort(:) = enswei(:)
            ENDIF

            ! Sort input ensemble
            ens_sort(:) = ens(ji,:)
            CALL heapsort2(ens_sort,enswei_sort)

            ! Compute cumulated weights
            DO jm=2,jpm
              enswei_ave(jm)=(enswei_sort(jm-1)+enswei_sort(jm))/2
            ENDDO

            enswei_sort(1) = 0.0
            DO jm=2,jpm
              enswei_sort(jm)=enswei_sort(jm-1)+enswei_ave(jm)
            ENDDO

            ! Rescale between 0 and 1
            enswei_sort(:) = enswei_sort(:) / enswei_sort(jpm)

            ! Get quantiles
            CALL get_quantiles( qua(ji,:), ens_sort, quadef, enswei_sort )

          ENDDO

          deallocate(enswei_sort,enswei_ave)

        ENDIF

        deallocate(ens_sort)

        END SUBROUTINE ens_quantiles_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ens_quantiles_variable( qua, ens, quadef, enswei )
!----------------------------------------------------------------------
!                  ***  ens_quantiles  ***
! 
! ** Purpose :   compute ensemble quantiles
!
! ** Arguments :
!         qua      : ensemble quantiles
!         ens      : input ensemble
!         quadef   : definition of the quantiles used for anamorphosis
!         enswei   : weight of ensemble members (optional)
!                    (default: all members have the same weight)
!         ensweiloc : local weight of ensemble members (optional)
!                     (default: all members have the same weight everywhere)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quadef 
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: enswei

        REAL(KIND=8), DIMENSION(:), allocatable :: ens_sort, enswei_sort, enswei_ave
        INTEGER :: jpm,jpq,jm,jq,allocstat

        jpm = SIZE(ens,1)
        jpq = SIZE(quadef,1)
        IF (SIZE(qua,1).NE.jpq) STOP 'Inconsistent size in ens_quantiles'
        IF (present(enswei)) THEN
          IF (SIZE(enswei,1).NE.jpm) STOP 'Inconsistent size in ens_quantiles'
        ENDIF

        IF (MINVAL(quadef).LT.0.0) STOP 'Bad quantile definition in ens_quantiles'
        IF (MAXVAL(quadef).GT.1.0) STOP 'Bad quantile definition in ens_quantiles'

        allocate ( ens_sort(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'

        IF (.NOT.present(enswei)) THEN

          ! Sort input ensemble
          ens_sort(:) = ens(:)
          CALL heapsort(ens_sort)

          ! Get quantiles
          CALL get_quantiles( qua, ens_sort, quadef )

        ELSE

          allocate ( enswei_sort(1:jpm), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'
          allocate ( enswei_ave(1:jpm), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in ens_quantiles'

          ! Sort input ensemble
          ens_sort(:) = ens(:)
          enswei_sort(:) = enswei(:)
          CALL heapsort2(ens_sort,enswei_sort)

          ! Compute cumulated weights
          DO jm=2,jpm
            enswei_ave(jm)=(enswei_sort(jm-1)+enswei_sort(jm))/2
          ENDDO

          enswei_sort(1) = 0.0
          DO jm=2,jpm
            enswei_sort(jm)=enswei_sort(jm-1)+enswei_ave(jm)
          ENDDO

          ! Rescale between 0 and 1
          enswei_sort(:) = enswei_sort(:) / enswei_sort(jpm)

          ! Get quantiles
          CALL get_quantiles( qua(:), ens_sort, quadef, enswei_sort )

          deallocate(enswei_sort,enswei_ave)

        ENDIF

        deallocate(ens_sort)

        END SUBROUTINE ens_quantiles_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE get_quantiles( qua, ens_sort, quadef, enswei )
!----------------------------------------------------------------------
! ** Purpose :   compute ensemble quantiles
!
! ** Arguments :
!         qua      : ensemble quantiles
!         ens      : input (sorted) ensemble
!         quadef   : definition of the quantiles
!         enswei   : cumulated weight of ensemble members (optional)
!                    (default: all members have the same weight)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: qua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ens_sort
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quadef
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: enswei

        INTEGER :: jpm,jpq,jm,jq,il,ih,im
        REAL(KIND=8) :: rm, a, b

        jpm = SIZE(ens_sort,1)
        jpq = SIZE(quadef,1)
        IF (SIZE(qua,1).NE.jpq) STOP 'Inconsistent size in get_quantiles'
        IF (present(enswei)) THEN
          IF (SIZE(enswei,1).NE.jpm) STOP 'Inconsistent size in get_quantiles'
        ENDIF

        IF (.NOT.present(enswei)) THEN

!         For each required quantile
          DO jq = 1,jpq

            rm = 1+quadef(jq)*(jpm-1)
            jm = INT(rm)

            b  = rm - jm
            a  = 1.0 - b

            qua(jq)=a*ens_sort(jm)
            IF (b.GT.0.) THEN
              qua(jq)=qua(jq)+b*ens_sort(jm+1)
            ENDIF

          ENDDO

        ELSE

!         For each required quantile
          DO jq = 1,jpq

!           Find rank by dichotomy in cumulated weights
            il = 1
            ih = jpm
            DO WHILE (il.LT.ih-1)
              im = (il+ih)/2
              IF (quadef(jq).GT.enswei(im)) THEN
                il = im
              ELSE
                ih = im
              ENDIF
            ENDDO
            jm = il
!           Compute interpolation weights (a,b)
            IF ( enswei(il+1) == enswei(il) ) THEN
              b = 0.
            ELSE
              b = ( quadef(jq) - enswei(il) ) / ( enswei(il+1) - enswei(il) )
            ENDIF
            a  = 1.0 - b
!           Interpolate to get the required quantile
            qua(jq)=a*ens_sort(jm)+b*ens_sort(jm+1)

          ENDDO

        ENDIF

        END SUBROUTINE get_quantiles
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Heapsort Fortran routines from Rosetta code:
! https://rosettacode.org/wiki/Sorting_algorithms/Heapsort
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
subroutine heapsort(a)
 
   real(kind=8), dimension(0:), intent(inout) :: a
   integer :: start, n, bottom
   real(kind=8) :: temp
 
   n = size(a)
   do start = (n - 2) / 2, 0, -1
     call siftdown(a, start, n);
   end do
 
   do bottom = n - 1, 1, -1
     temp = a(0)
     a(0) = a(bottom)
     a(bottom) = temp;
     call siftdown(a, 0, bottom)
   end do
 
end subroutine heapsort
 
subroutine siftdown(a, start, bottom)
 
  real(kind=8), intent(in out) :: a(0:)
  integer, intent(in) :: start, bottom
  integer :: child, root
  real(kind=8) :: temp
 
  root = start
  do while(root*2 + 1 < bottom)
    child = root * 2 + 1
 
    if (child + 1 < bottom) then
      if (a(child) < a(child+1)) child = child + 1
    end if
 
    if (a(root) < a(child)) then
      temp = a(child)
      a(child) = a (root)
      a(root) = temp
      root = child
    else
      return
    end if  
  end do      
 
end subroutine siftdown

subroutine heapsort2(a,b)

   real(kind=8), dimension(0:), intent(in out) :: a
   real(kind=8), dimension(0:), intent(in out) :: b
   integer :: start, n, bottom
   real(kind=8) :: temp

   n = size(a)
   do start = (n - 2) / 2, 0, -1
     call siftdown2(a, b, start, n);
   end do

   do bottom = n - 1, 1, -1
     temp = a(0)
     a(0) = a(bottom)
     a(bottom) = temp;
     temp = b(0)
     b(0) = b(bottom)
     b(bottom) = temp;
     call siftdown2(a, b, 0, bottom)
   end do

end subroutine heapsort2

subroutine siftdown2(a, b, start, bottom)

  real(kind=8), intent(in out) :: a(0:)
  real(kind=8), intent(in out) :: b(0:)
  integer, intent(in) :: start, bottom
  integer :: child, root
  real(kind=8) :: temp

  root = start
  do while(root*2 + 1 < bottom)
    child = root * 2 + 1

    if (child + 1 < bottom) then
      if (a(child) < a(child+1)) child = child + 1
    end if

    if (a(root) < a(child)) then
      temp = a(child)
      a(child) = a (root)
      a(root) = temp
      temp = b(child)
      b(child) = b (root)
      b(root) = temp
      root = child
    else
      return
    end if
  end do

end subroutine siftdown2
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Heapsort Fortran routines from Rosetta code:
! https://rosettacode.org/wiki/Sorting_algorithms/Heapsort
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ensdam_anaqua
