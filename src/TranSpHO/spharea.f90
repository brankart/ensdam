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
!                        MODULE SPHAREA
!
!---------------------------------------------------------------------
! Computation of spherical areas
! by Jean-Michel Brankart, December 2017
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
MODULE ensdam_spharea
        IMPLICIT NONE
        PRIVATE

        PUBLIC mesh_area

        REAL(KIND=8), PARAMETER :: twopi=2*3.1415926535897932384626
        REAL(KIND=8), PARAMETER :: deg2rad=twopi/360.

CONTAINS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  -----------------------------------------------------------------
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE mesh_area (lon,lat,area)
!---------------------------------------------------------------------
!
!  Purpose : Compute sphercial surface element for each grid point
!
!  Method : Compute surface of each quadrangle using spherical trigonometry
!
!  Input : lon,lat : latitude and longitude of each grid point
!  Output : area : surface associated to each grid point (sphere = 1)
!
!---------------------------------------------------------------------
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: lon, lat
          REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: area

          INTEGER :: jpi,jpj,ji,jj,jpigrid,jpjgrid
          REAL(KIND=8) :: lat1,lat2,lat3,lat4
          REAL(KIND=8) :: lon1,lon2,lon3,lon4
          REAL(KIND=8) :: c12, c23, c34, c41, c13, s, e, a123, a143, a

          jpi=SIZE(area,1)
          jpj=SIZE(area,2)
          jpigrid=SIZE(lon,1)
          jpjgrid=SIZE(lon,2)
          IF (jpigrid.NE.SIZE(lat,1)) STOP 'Bad grid size in spharea'
          IF (jpjgrid.NE.SIZE(lat,2)) STOP 'Bad grid size in spharea'
          IF (jpigrid.NE.SIZE(area,1)+1) STOP 'Bad grid size in spharea'
          IF (jpjgrid.NE.SIZE(area,2)+1) STOP 'Bad grid size in spharea'

          DO ji=1,jpi
          DO jj=1,jpj
            ! coordinates of the vertices of the quadrangle
            lon1 = lon(ji,jj) * deg2rad
            lon2 = lon(ji+1,jj) * deg2rad
            lon3 = lon(ji+1,jj+1) * deg2rad
            lon4 = lon(ji,jj+1) * deg2rad
            lat1 = ( 90. - lat(ji,jj) ) * deg2rad
            lat2 = ( 90. - lat(ji+1,jj) ) * deg2rad
            lat3 = ( 90. - lat(ji+1,jj+1) ) * deg2rad
            lat4 = ( 90. - lat(ji,jj+1) ) * deg2rad
            ! distance between vertices
            c12 = ACOS ( MIN( SIN(lat1)*SIN(lat2)*COS(lon1-lon2) + COS(lat1)*COS(lat2) , 1.) )
            c23 = ACOS ( MIN( SIN(lat2)*SIN(lat3)*COS(lon2-lon3) + COS(lat2)*COS(lat3) , 1.) )
            c34 = ACOS ( MIN( SIN(lat3)*SIN(lat4)*COS(lon3-lon4) + COS(lat3)*COS(lat4) , 1.) )
            c41 = ACOS ( MIN( SIN(lat4)*SIN(lat1)*COS(lon4-lon1) + COS(lat4)*COS(lat1) , 1.) )
            c13 = ACOS ( MIN( SIN(lat1)*SIN(lat3)*COS(lon1-lon3) + COS(lat1)*COS(lat3) , 1.) )
            ! area of 1-2-3 triangle
            s = 0.5 * ( c12 + c23 + c13 )
            e = sqrt(MAX(tan(s/2)*tan((s-c12)/2)*tan((s-c23)/2)*tan((s-c13)/2),0.))
            a123 = 4.0 * ATAN(e)
            ! area of 1-4-3 triangle
            s = 0.5 * ( c34 + c41 + c13 )
            e = sqrt(MAX(tan(s/2)*tan((s-c34)/2)*tan((s-c41)/2)*tan((s-c13)/2),0.))
            a143 = 4.0 * ATAN(e)
            ! total area of quadrangle (sphere = 4pi)
            a = a123 + a143

            area(ji,jj) = 0.5 * a / twopi
          ENDDO
          ENDDO

        END SUBROUTINE mesh_area
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_spharea
