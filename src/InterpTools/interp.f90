! Copyright: CNRS - Université de Grenoble

! Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
!                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio

! Jean-Michel.Brankart@hmg.inpg.fr

! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".

!-----------------------------------------------------------------------

!                        MODULE INTERP

!-----------------------------------------------------------------------
! Interpolation in gridded data
!   - 1d or 2D input grid
!   - cartesian or spherical coordinates

! by Jean-Michel Brankart, June 2018
! ----------------------------------------------------------------------
      MODULE ensdam_interp
        IMPLICIT NONE
        PRIVATE

        PUBLIC grid1D_locate
        PUBLIC grid1D_interp
        PUBLIC grid2D_init
        PUBLIC grid2D_locate
        PUBLIC grid2D_interp
        PUBLIC sph_dist

        LOGICAL, SAVE :: periodic = .TRUE. ! periodic search in 2D grid
        CHARACTER(len=20) :: grid_type = 'cartesian' ! type of 2D grid

        INTEGER, SAVE :: jpi, jpj ! 2D grid size
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: xg, yg ! 2D grid arrays
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: viewed ! grid points already checked
        INTEGER, SAVE :: jiloc, jjloc ! Current grid cell
        REAL(kind=8), SAVE :: dcellmin,dcellold ! distance of cells to observation
        LOGICAL, SAVE :: located,stopped ! observation is located or not
        INTEGER, SAVE :: ngrow ! how many times is distance to observation is growing

        REAL(KIND=8), PARAMETER :: twopi=2*3.1415926535897932384626
        REAL(KIND=8), PARAMETER :: deg2rad=twopi/360.

      CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE grid2D_init( kxg, kyg, gtype )
!----------------------------------------------------------------------
!                  ***  grid2D_init  ***
! 
! ** Purpose :   define the grid mesh
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: kxg, kyg
        CHARACTER(len=*), INTENT( in ) :: gtype

        INTEGER :: allocok

! --- Define grid size
        jpi=SIZE(kxg,1)
        jpj=SIZE(kxg,2)

        grid_type = gtype

        IF ( (SIZE(kyg,1).NE.jpi).OR.(SIZE(kyg,2).NE.jpj) ) THEN
          PRINT *, 'Incorrect grid size in grid2D_init'
        ENDIF

! --- Allocate grid arrays
        allocate(xg(jpi,jpj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in grid2D_init'
        allocate(yg(jpi,jpj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in grid2D_init'
        allocate(viewed(jpi,jpj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in grid2D_init'

! --- initialize grid arrays
        xg(:,:) = kxg(:,:)
        yg(:,:) = kyg(:,:)
        viewed (:,:) = .FALSE.

! --- set starting research cell for 1st observation
        jiloc = jpi/2
        jjloc = jpj/2

        END SUBROUTINE grid2D_init
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION grid1D_locate( kgrid, kx, ki )
!----------------------------------------------------------------------
!                  ***  grid1D_locate ***
! 
! ** Purpose :   locate data point in grid mesh
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kgrid
        REAL(KIND=8), INTENT( in ) :: kx
        INTEGER, INTENT( out ) :: ki
        LOGICAL :: grid1D_locate

        INTEGER :: jpi,ibelow,imid,iabove

        jpi=SIZE(kgrid,1)

        located = ( kx.GE.kgrid(1) ) .AND. ( kx.LE.kgrid(jpi) )

        IF (located) THEN
          ibelow=1
          imid=jpi/2
          iabove=jpi

          DO WHILE (ibelow.NE.imid)
            IF (kx.GE.kgrid(imid)) THEN
              ibelow=imid
            ELSE
              iabove=imid
            ENDIF
            imid=ibelow+(iabove-ibelow)/2
          ENDDO

          ki=imid
        ENDIF

        grid1D_locate = located

        END FUNCTION grid1D_locate
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION grid2D_locate( kx, ky, ki, kj )
!----------------------------------------------------------------------
!                  ***  grid2D_locate ***
! 
! ** Purpose :   locate data point in grid mesh
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: kx, ky
        INTEGER, INTENT( out ) :: ki, kj
        LOGICAL :: grid2D_locate

        viewed(:,:) = .FALSE.
        stopped  = .FALSE.
        ngrow = 0
        dcellold = huge(dcellold)

        located = incell(jiloc,jjloc,kx,ky)
        viewed(jiloc,jjloc) = .TRUE.

! --- follow a path accross the mesh towards the observation cell
        DO WHILE ( (.NOT.located) .AND. (.NOT.stopped) )

! --- move to neighbour cell nearest to obseration
          CALL move_to_neighbour_cell( kx, ky )

! --- test new cell
          located = incell(jiloc,jjloc,kx,ky)
          viewed(jiloc,jjloc) = .TRUE.

        ENDDO

        IF (located) THEN
          ki = jiloc
          kj = jjloc
        ENDIF

        grid2D_locate = located

        END FUNCTION grid2D_locate
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE grid1D_interp( kgrid, kx, ki, w )
!----------------------------------------------------------------------
!                  ***  grid1D_interp  ***
! 
! ** Purpose :   compute interpolation weights
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: kgrid
        REAL(KIND=8), INTENT( in ) :: kx
        INTEGER, INTENT( in ) :: ki
        REAL(KIND=8), INTENT( out ) :: w

        w=(kx-kgrid(ki))/(kgrid(ki+1)-kgrid(ki))

        END SUBROUTINE grid1D_interp
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE grid2D_interp( kx, ky, ki, kj, kw )
!----------------------------------------------------------------------
!                  ***  grid2D_interp  ***
! 
! ** Purpose :   compute interpolation weights
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: kx, ky
        INTEGER, INTENT( in ) :: ki, kj
        REAL(KIND=8), DIMENSION(2,2), INTENT( out ) :: kw

        REAL(KIND=8) :: r,s
        REAL(KIND=8), dimension(4) :: jix,jjy

        REAL(KIND=8) :: lat1,lat2,lat3,lat4,olat
        REAL(KIND=8) :: lon1,lon2,lon3,lon4,olon
        REAL(KIND=8) :: a12,a23,a34,a41,atot
        REAL(KIND=8) :: d1, d2, d3, d4
        REAL(KIND=8) :: kw1, kw2

        SELECT CASE(grid_type)
        CASE ('cartesian')

! --- position of the cell points
! 4---------3
! |   obs   |
! 1---------2
          jix (1) = xg(ki,kj)
          jix (2) = xg(ki+1,kj)
          jix (3) = xg(ki+1,kj+1)
          jix (4) = xg(ki,kj+1)
          jjy (1) = yg(ki,kj)
          jjy (2) = yg(ki+1,kj)
          jjy (3) = yg(ki+1,kj+1)
          jjy (4) = yg(ki,kj+1)

! --- remap cell on a unit square
          CALL mkxytors(kx,ky,jix,jjy,r,s)

! --- compute interpolation weights
          kw(1,1) = 0.25 * ( 1 - r ) * ( 1 - s )
          kw(2,1) = 0.25 * ( 1 + r ) * ( 1 - s )
          kw(2,2) = 0.25 * ( 1 + r ) * ( 1 + s )
          kw(1,2) = 0.25 * ( 1 - r ) * ( 1 + s )

        CASE('spherical')

! --- position of the cell points
! 4---------3
! |   obs   |
! 1---------2
          lon1 = xg(ki,kj) * deg2rad
          lon2 = xg(ki+1,kj) * deg2rad
          lon3 = xg(ki+1,kj+1) * deg2rad
          lon4 = xg(ki,kj+1) * deg2rad
          lat1 = ( 90. - yg(ki,kj) ) * deg2rad
          lat2 = ( 90. - yg(ki+1,kj) ) * deg2rad
          lat3 = ( 90. - yg(ki+1,kj+1) ) * deg2rad
          lat4 = ( 90. - yg(ki,kj+1) ) * deg2rad

! --- position of observation point
          olon = kx * deg2rad
          olat = ( 90. - ky ) * deg2rad

! --- distance to the cell vertices: o1, o2, o3, o4
          d1 = sph_dist(olon,lon1,olat,lat1)
          d2 = sph_dist(olon,lon2,olat,lat2)
          d3 = sph_dist(olon,lon3,olat,lat3)
          d4 = sph_dist(olon,lon4,olat,lat4)

! --- area of triangles: o12, o23, o34, o41
          a12 = sph_area(olon,lon1,lon2,olat,lat1,lat2)
          a23 = sph_area(olon,lon2,lon3,olat,lat2,lat3)
          a34 = sph_area(olon,lon3,lon4,olat,lat3,lat4)
          a41 = sph_area(olon,lon4,lon1,olat,lat4,lat1)

          kw(1,1) = a23*a34*a41 * d2/(d1+d2) + a12*a23*a34 * d4/(d1+d4)
          kw(2,1) = a23*a34*a41 * d1/(d1+d2) + a12*a34*a41 * d3/(d2+d3)
          kw(2,2) = a12*a34*a41 * d2/(d2+d3) + a12*a23*a41 * d4/(d4+d3)
          kw(1,2) = a12*a23*a41 * d3/(d4+d3) + a12*a23*a34 * d1/(d1+d4)

          atot = a23*a34*a41 + a12*a34*a41 + a12*a23*a41 + a12*a23*a34
          kw(:,:) = kw(:,:) / atot

        CASE('spherical_old')

! --- position of the cell points
! 4---------3
! |   obs   |
! 1---------2
          lon1 = xg(ki,kj) * deg2rad
          lon2 = xg(ki+1,kj) * deg2rad
          lon3 = xg(ki+1,kj+1) * deg2rad
          lon4 = xg(ki,kj+1) * deg2rad
          lat1 = ( 90. - yg(ki,kj) ) * deg2rad
          lat2 = ( 90. - yg(ki+1,kj) ) * deg2rad
          lat3 = ( 90. - yg(ki+1,kj+1) ) * deg2rad
          lat4 = ( 90. - yg(ki,kj+1) ) * deg2rad

! --- position of observation point
          olon = kx * deg2rad
          olat = ( 90. - ky ) * deg2rad

! --- area of triangles o12, o23, o34, o41
          a12 = sph_area(olon,lon1,lon2,olat,lat1,lat2)
          a23 = sph_area(olon,lon2,lon3,olat,lat2,lat3)
          a34 = sph_area(olon,lon3,lon4,olat,lat3,lat4)
          a41 = sph_area(olon,lon4,lon1,olat,lat4,lat1)
          atot = ( a12 + a23 + a34 + a41 ) * 2.0

! --- compute interpolation weights
          kw1 = a23
          atot = a41 + a23
          IF (atot.GT.0.0) kw1 = kw1 / atot

          kw2 = a34
          atot = a12 + a34
          IF (atot.GT.0.0) kw2 = kw2 / atot

          kw(1,1) = kw1 * kw2
          kw(2,1) = ( 1.0 - kw1 ) * kw2
          kw(2,2) = ( 1.0 - kw1 ) * ( 1.0 - kw2 )
          kw(1,2) = kw1 * ( 1.0 - kw2 )

        CASE DEFAULT
          STOP 'invalid grid type in mod_interp'
        END SELECT

        END SUBROUTINE grid2D_interp
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE move_to_neighbour_cell( kx, ky )
!----------------------------------------------------------------------
!                  ***  move_to_neighbour_cell  ***
! 
! ** Purpose :   move to neighbour cell nearest to observation
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: kx, ky

        INTEGER, dimension(4) :: ivic,jvic
        INTEGER :: icell
        REAL(kind=8) :: dcell

! --- define neighbour cells
        IF (periodic) THEN 
          ivic(1) = jiloc-1
          IF (jiloc-1.EQ.0) ivic(1) = jpi-1
          jvic(1) = jjloc
          ivic(2) = jiloc
          jvic(2) = jjloc-1
          IF (jjloc-1.EQ.0) jvic(2) = jpj-1
          ivic(3) = jiloc+1
          IF (jiloc+1.EQ.jpi) ivic(3) = 1
          jvic(3) = jjloc
          ivic(4) = jiloc
          jvic(4) = jjloc+1
          IF (jjloc+1.EQ.jpj) jvic(4) = 1
        ELSE
          ivic(1) = MAX(1,jiloc-1)
          jvic(1) = jjloc
          ivic(2) = jiloc
          jvic(2) = MAX(1,jjloc-1)
          ivic(3) = MIN(jiloc+1,jpi-1)
          jvic(3) = jjloc
          ivic(4) = jiloc
          jvic(4) = MIN(jjloc+1,jpj-1)
        ENDIF

! --- find the one that is closest to observation
        dcellmin = huge(dcellmin)
        DO icell = 1,4
          IF ( .NOT.viewed(ivic(icell),jvic(icell)) ) THEN
            dcell = cell_distance(ivic(icell),jvic(icell),kx,ky)
            IF (dcell.LT.dcellmin) THEN
              dcellmin = dcell
              jiloc = ivic(icell)
              jjloc = jvic(icell)
            ENDIF
          ENDIF
        ENDDO

! --- detect abnormal end of path
        stopped = dcellmin.EQ.huge(dcellmin) ! no more possible move
        stopped = stopped .OR. (ngrow.GT.3)  ! path going away

! --- store information when path is going away from observation
        IF (dcellmin.GT.dcellold) ngrow = ngrow + 1
        IF (dcellmin.LT.dcellold) ngrow = 0
        dcellold = dcellmin

        END SUBROUTINE move_to_neighbour_cell
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION cell_distance(ki,kj,kx,ky)
!----------------------------------------------------------------------
!                  ***  cell_distance  ***
! 
! ** Purpose :   compute distance between grid cell and observation
!----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT( in ) :: ki, kj
        REAL(KIND=8), INTENT( in ) :: kx, ky
        REAL(KIND=8) :: cell_distance

        REAL(KIND=8) :: xbaryc, ybaryc, dx, dy, d

        REAL(KIND=8) :: lat1,lat2,lat3,lat4,olat
        REAL(KIND=8) :: lon1,lon2,lon3,lon4,olon

        SELECT CASE(grid_type)
        CASE ('cartesian')

          xbaryc = xg(ki,kj) +  xg(ki+1,kj) &
     &           + xg(ki,kj+1) + xg(ki+1,kj+1)
          ybaryc = yg(ki,kj) +  yg(ki+1,kj) &
     &           + yg(ki,kj+1) + yg(ki+1,kj+1)

          xbaryc = xbaryc * 0.25 ; ybaryc = ybaryc * 0.25

          dx = kx - xbaryc ; dy = ky - ybaryc

          cell_distance = dx * dx + dy * dy

        CASE('spherical')

          lon1 = xg(ki,kj) * deg2rad
          lon2 = xg(ki+1,kj) * deg2rad
          lon3 = xg(ki+1,kj+1) * deg2rad
          lon4 = xg(ki,kj+1) * deg2rad
          lat1 = ( 90. - yg(ki,kj) ) * deg2rad
          lat2 = ( 90. - yg(ki+1,kj) ) * deg2rad
          lat3 = ( 90. - yg(ki+1,kj+1) ) * deg2rad
          lat4 = ( 90. - yg(ki,kj+1) ) * deg2rad

          olon = kx * deg2rad
          olat = ( 90. - ky ) * deg2rad

          cell_distance = HUGE(cell_distance)
          d = sph_dist(olon,lon1,olat,lat1)
          cell_distance = MIN(cell_distance,d)
          d = sph_dist(olon,lon2,olat,lat2)
          cell_distance = MIN(cell_distance,d)
          d = sph_dist(olon,lon3,olat,lat3)
          cell_distance = MIN(cell_distance,d)
          d = sph_dist(olon,lon4,olat,lat4)
          cell_distance = MIN(cell_distance,d)

        CASE DEFAULT
          STOP 'invalid grid type in mod_interp'
        END SELECT

        END FUNCTION cell_distance
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION incell(ki,kj,kx,ky)
!----------------------------------------------------------------------
!                  ***  incell  ***
! 
! ** Purpose :   check if observation point is inside current cell
!----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT( in ) :: ki, kj
        REAL(KIND=8), INTENT( in ) :: kx, ky
        LOGICAL :: incell

        REAL(KIND=8) :: lat1,lat2,lat3,lat4,olat
        REAL(KIND=8) :: lon1,lon2,lon3,lon4,olon
        REAL(KIND=8) :: atri12,atri23,atri34,atri41
        REAL(KIND=8) :: atri13,atri,atot

        SELECT CASE(grid_type)
        CASE ('cartesian')

          incell = &
     &    sameside(kx, xg(ki+1,kj+1), xg(ki,kj),   xg(ki+1,kj),&
     &             ky, yg(ki+1,kj+1), yg(ki,kj),   yg(ki+1,kj) )&
     &    .AND.&
     &    sameside(kx, xg(ki,kj+1),   xg(ki+1,kj), xg(ki+1,kj+1),&
     &             ky, yg(ki,kj+1),   yg(ki+1,kj), yg(ki+1,kj+1) )&
     &    .AND.&
     &    sameside(kx, xg(ki,kj),   xg(ki+1,kj+1), xg(ki,kj+1),&
     &             ky, yg(ki,kj),   yg(ki+1,kj+1), yg(ki,kj+1) )&
     &    .AND.&
     &    sameside(kx, xg(ki+1,kj),   xg(ki,kj+1), xg(ki,kj),&
     &             ky, yg(ki+1,kj),   yg(ki,kj+1), yg(ki,kj) )

        CASE('spherical')

! --- position of the cell points
! 4---------3
! |   obs   |
! 1---------2
          lon1 = xg(ki,kj) * deg2rad
          lon2 = xg(ki+1,kj) * deg2rad
          lon3 = xg(ki+1,kj+1) * deg2rad
          lon4 = xg(ki,kj+1) * deg2rad
          lat1 = ( 90. - yg(ki,kj) ) * deg2rad
          lat2 = ( 90. - yg(ki+1,kj) ) * deg2rad
          lat3 = ( 90. - yg(ki+1,kj+1) ) * deg2rad
          lat4 = ( 90. - yg(ki,kj+1) ) * deg2rad

! --- position of observation point
          olon = kx * deg2rad
          olat = ( 90. - ky ) * deg2rad

! --- area of triangles 123, o12, o23, o13
          atri   = sph_area(lon1,lon2,lon3,lat1,lat2,lat3)
          atri12 = sph_area(olon,lon1,lon2,olat,lat1,lat2)
          atri23 = sph_area(olon,lon2,lon3,olat,lat2,lat3)
          atri13 = sph_area(olon,lon4,lon1,olat,lat4,lat1)
          incell = atri12+atri23+atri13 .LE. atri

! --- area of triangles 134, o34, o41
          IF (.NOT.incell ) THEN
            atri   = sph_area(lon1,lon3,lon4,lat1,lat3,lat4)
            atri34 = sph_area(olon,lon3,lon4,olat,lat3,lat4)
            atri41 = sph_area(olon,lon4,lon1,olat,lat4,lat1)
            incell = atri34+atri41+atri13 .LE. atri
          ENDIF

        CASE DEFAULT
          STOP 'invalid grid type in mod_interp'
        END SELECT

        END FUNCTION incell
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION sph_dist(lon1,lon2,lat1,lat2)
!----------------------------------------------------------------------
!                  ***  sph_dist  ***
! 
! ** Purpose : compute distance along the sphere
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: lon1,lon2,lat1,lat2
        REAL(KIND=8) :: sph_dist

        sph_dist = ACOS ( MIN( SIN(lat1)*SIN(lat2)*COS(lon1-lon2) &
     &                    + COS(lat1)*COS(lat2) , 1.) )

        END FUNCTION sph_dist
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION sph_area(lon1,lon2,lon3,lat1,lat2,lat3)
!----------------------------------------------------------------------
!                  ***  sph_area  ***
! 
! ** Purpose : compute surface of spherical triangle
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: lon1,lon2,lon3,lat1,lat2,lat3
        REAL(KIND=8) :: sph_area

        REAL(KIND=8) :: c12, c23, c31, s, e

! --- distance between vertices
        c12 = sph_dist(lon1,lon2,lat1,lat2)
        c23 = sph_dist(lon2,lon3,lat2,lat3)
        c31 = sph_dist(lon3,lon1,lat3,lat1)

! --- area of triangle
        s = 0.5 * ( c12 + c23 + c31 )
        e = tan(s/2)*tan((s-c12)/2)*tan((s-c23)/2)*tan((s-c31)/2)
        e = SQRT ( MAX(e,0.) )
        sph_area = 4.0 * ATAN(e)

        END FUNCTION sph_area
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION sameside(fxa,fxb,fx1,fx2,fya,fyb,fy1,fy2)
!----------------------------------------------------------------------
!                  ***  sameside  ***
! 
! ** Purpose :   check if 2 points are on the same side of a straight line
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: fxa,fxb,fx1,fx2,fya,fyb,fy1,fy2
        LOGICAL :: sameside

        REAL(KIND=8) :: fsidea, fsideb

        fsidea = (fya-fy1)*(fx2-fx1)-(fxa-fx1)*(fy2-fy1)
        fsideb = (fyb-fy1)*(fx2-fx1)-(fxb-fx1)*(fy2-fy1)
        sameside = ( fsidea * fsideb .GE. 0. )

        END FUNCTION sameside
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE mkxytors(jox,joy,jix,jjy,r,s)

!     Inverse the nonlinear quadrangular mapping 
!     using a bidimensional Newton-Raphson method

        IMPLICIT NONE

        REAL(KIND=8), intent(in) :: jox,joy
        REAL(KIND=8), dimension(4), intent(in) :: jix,jjy
        REAL(KIND=8), intent(out) :: r,s

        REAL(KIND=8) :: rho,epsilon,rr,ss
        REAL(KIND=8), dimension(2,2) :: aa,aai,aaj
        REAL(KIND=8), dimension(2) :: bb,tt,cc

        INTEGER :: iter

        epsilon = 0.0001

        aa(1,1) = 0.25 * ( jix(2) + jix(3) - jix(1) - jix(4) )
        aa(1,2) = 0.25 * ( jix(3) + jix(4) - jix(1) - jix(2) )
        aa(2,1) = 0.25 * ( jjy(2) + jjy(3) - jjy(1) - jjy(4) )
        aa(2,2) = 0.25 * ( jjy(3) + jjy(4) - jjy(1) - jjy(2) )

        bb(1) = jox - 0.25 * ( jix(1) + jix(2) + jix(3) + jix(4) )
        bb(2) = joy - 0.25 * ( jjy(1) + jjy(2) + jjy(3) + jjy(4) )

        tt(1) = 0.25 * ( jix(1) + jix(3) - jix(2) - jix(4) )
        tt(2) = 0.25 * ( jjy(1) + jjy(3) - jjy(2) - jjy(4) )

        rr = 0.
        ss = 0.
        LOOP1 : DO iter = 1,100
           aaj(1,1) = aa(1,1) + ss * tt(1)
           aaj(1,2) = aa(1,2) + rr * tt(1)
           aaj(2,1) = aa(2,1) + ss * tt(2)
           aaj(2,2) = aa(2,2) + rr * tt(2)

           rho = aaj(1,1) * aaj(2,2) - aaj(1,2) * aaj(2,1)

           aai(1,1) = aaj(2,2) / rho
           aai(1,2) = - aaj(1,2) / rho
           aai(2,1) = - aaj(2,1) / rho
           aai(2,2) = aaj(1,1) / rho

           cc(1) = rr * ss * tt(1) - bb(1)
           cc(2) = rr * ss * tt(2) - bb(2)
           cc(1) = cc(1) + aa(1,1) * rr + aa(1,2) * ss
           cc(2) = cc(2) + aa(2,1) * rr + aa(2,2) * ss

           r = rr - aai(1,1) * cc(1) - aai(1,2) * cc(2)
           s = ss - aai(2,1) * cc(1) - aai(2,2) * cc(2)

           IF ( (ABS(r-rr).LT.epsilon) .AND. (ABS(s-ss).LT.epsilon) ) THEN
              exit LOOP1
           ELSE
              rr = r
              ss = s
           ENDIF
        ENDDO LOOP1

        IF (r.GT.1.) r = 1.
        IF (s.GT.1.) s = 1.
        IF (r.LT.-1.) r = -1.
        IF (s.LT.-1.) s = -1.

        END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE mkrstoxy(jox,joy,jix,jjy,r,s)

!     Quadrangular mapping 

        IMPLICIT NONE
        REAL(KIND=8), intent(out) :: jox,joy
        REAL(KIND=8), dimension(4), intent(in) :: jix,jjy
        REAL(KIND=8), intent(in) :: r,s

        REAL(KIND=8) :: fn1,fn2,fn3,fn4

        fn1 = 0.25 * ( 1 - r ) * ( 1 - s )
        fn2 = 0.25 * ( 1 + r ) * ( 1 - s )
        fn3 = 0.25 * ( 1 + r ) * ( 1 + s )
        fn4 = 0.25 * ( 1 - r ) * ( 1 + s )

        jox = fn1 * jix(1) + fn2 * jix(2) + fn3 * jix(3) + fn4 * jix(4)
        joy = fn1 * jjy(1) + fn2 * jjy(2) + fn3 * jjy(3) + fn4 * jjy(4)

        END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE ensdam_interp
