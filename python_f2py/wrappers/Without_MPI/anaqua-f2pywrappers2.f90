!     -*- f90 -*-
!     This file is autogenerated with f2py (version:1.21.6)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_ensdam_anaqua_ens_quantiles_vector (qua, ens, &
     &quadef, enswei, ensweiloc, f2py_ens_d0, f2py_ens_d1, f2py_quadef_d&
     &0)
      use ensdam_anaqua, only : ens_quantiles_vector
      integer f2py_ens_d0
      integer f2py_ens_d1
      integer f2py_quadef_d0
      real(kind=8) ens(f2py_ens_d0,f2py_ens_d1)
      real(kind=8) quadef(f2py_quadef_d0)
      real(kind=8) qua(size(ens,1),size(quadef,1))
      real(kind=8) enswei(size(ens,2))
      real(kind=8) ensweiloc(size(ens,1),size(ens,2))
      call ens_quantiles_vector(qua, ens, quadef, enswei, ensweiloc)
      end subroutine f2pywrap_ensdam_anaqua_ens_quantiles_vector
      subroutine f2pywrap_ensdam_anaqua_ens_quantiles_variable (qua, ens&
     &, quadef, enswei, f2py_qua_d0, f2py_ens_d0)
      use ensdam_anaqua, only : ens_quantiles_variable
      integer f2py_qua_d0
      integer f2py_ens_d0
      real(kind=8) qua(f2py_qua_d0)
      real(kind=8) ens(f2py_ens_d0)
      real(kind=8) quadef(size(qua,1))
      real(kind=8) enswei(size(ens,1))
      call ens_quantiles_variable(qua, ens, quadef, enswei)
      end subroutine f2pywrap_ensdam_anaqua_ens_quantiles_variable
      
      subroutine f2pyinitensdam_anaqua(f2pysetupfunc)
      use ensdam_anaqua, only : heapsort
      use ensdam_anaqua, only : heapsort2
      interface 
      subroutine f2pywrap_ensdam_anaqua_ens_quantiles_vector (qua, ens, &
     &quadef, enswei, ensweiloc, f2py_ens_d0, f2py_ens_d1, f2py_quadef_d&
     &0)
      integer f2py_ens_d0
      integer f2py_ens_d1
      integer f2py_quadef_d0
      real(kind=8) ens(f2py_ens_d0,f2py_ens_d1)
      real(kind=8) quadef(f2py_quadef_d0)
      real(kind=8) qua(size(ens,1),size(quadef,1))
      real(kind=8) enswei(size(ens,2))
      real(kind=8) ensweiloc(size(ens,1),size(ens,2))
      end subroutine f2pywrap_ensdam_anaqua_ens_quantiles_vector 
      subroutine f2pywrap_ensdam_anaqua_ens_quantiles_variable (qua, ens&
     &, quadef, enswei, f2py_qua_d0, f2py_ens_d0)
      integer f2py_qua_d0
      integer f2py_ens_d0
      real(kind=8) qua(f2py_qua_d0)
      real(kind=8) ens(f2py_ens_d0)
      real(kind=8) quadef(size(qua,1))
      real(kind=8) enswei(size(ens,1))
      end subroutine f2pywrap_ensdam_anaqua_ens_quantiles_variable
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_ensdam_anaqua_ens_quantiles_vector,f2p&
     &ywrap_ensdam_anaqua_ens_quantiles_variable,heapsort,heapsort2)
      end subroutine f2pyinitensdam_anaqua

