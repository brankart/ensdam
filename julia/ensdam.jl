#---------------------------------------------------------------------
# Copyright: CNRS - Université de Grenoble Alpes
#
# Contributors : Jean-Michel Brankart
#
# Jean-Michel.Brankart@univ-grenoble-alpes.fr
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
#---------------------------------------------------------------------

module ensdam

const libdir = "../lib"
const lib = joinpath(libdir,"libensdam") 

####################
# Wrapper to stoutil
####################

function cdf_gaussian(x)
  y = ccall((:ensdam_stoutil_mp_cdf_gaussian_,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function pdf_gaussian(x)
  y = ccall((:ensdam_stoutil_mp_pdf_gaussian_,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function logpdf_gaussian(x)
  y = ccall((:ensdam_stoutil_mp_logpdf_gaussian_,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function invcdf_gaussian(x)
  y = ccall((:ensdam_stoutil_mp_invcdf_gaussian_,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function cdf_gamma(a,x)
  y = ccall((:ensdam_stoutil_mp_cdf_gamma_,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function pdf_gamma(a,x)
  y = ccall((:ensdam_stoutil_mp_pdf_gamma_,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function logpdf_gamma(a,x)
  y = ccall((:ensdam_stoutil_mp_logpdf_gamma_,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function invcdf_gamma(a,x)
  y = ccall((:ensdam_stoutil_mp_invcdf_gamma_,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function cdf_beta(a,b,x)
  y = ccall((:ensdam_stoutil_mp_cdf_beta_,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function pdf_beta(a,b,x)
  y = ccall((:ensdam_stoutil_mp_pdf_beta_,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function logpdf_beta(a,b,x)
  y = ccall((:ensdam_stoutil_mp_logpdf_beta_,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function invcdf_beta(a,b,x)
  y = ccall((:ensdam_stoutil_mp_invcdf_beta_,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

end # end of module
