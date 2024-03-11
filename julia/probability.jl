module probability

const libdir = dirname(pathof(@__MODULE__))
const lib = joinpath(libdir,"ensdam") 

###################
# Wrapper to storng
###################

function kiss_uniform()
  x=Float64[1]
  ccall((:c_kiss_uniform,lib), Cvoid, (Ptr{Float64},), x)
  return x[1]
end

function kiss_gaussian()
  x=Float64[1]
  ccall((:c_kiss_gaussian,lib), Cvoid, (Ptr{Float64},), x)
  return x[1]
end

function kiss_gamma(k)
  k=Float64[k]
  x=Float64[1]
  ccall((:c_kiss_gamma,lib), Cvoid, (Ptr{Float64},Ptr{Float64}), x, k)
  return x[1]
end

function kiss_beta(a,b)
  a=Float64[a]
  b=Float64[b]
  x=Float64[1]
  ccall((:c_kiss_beta,lib), Cvoid, (Ptr{Float64},Ptr{Float64},Ptr{Float64}), x, a, b)
  return x[1]
end

####################
# Wrapper to stoutil
####################

function cdf_gaussian(x)
  y = ccall((:c_cdf_gaussian,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function pdf_gaussian(x)
  y = ccall((:c_pdf_gaussian,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function logpdf_gaussian(x)
  y = ccall((:c_logpdf_gaussian,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function invcdf_gaussian(x)
  y = ccall((:c_invcdf_gaussian,lib), Float64, (Ref{Float64},), Float64(x))
  return y
end

function cdf_gamma(a,x)
  y = ccall((:c_cdf_gamma,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function pdf_gamma(a,x)
  y = ccall((:c_pdf_gamma,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function logpdf_gamma(a,x)
  y = ccall((:c_logpdf_gamma,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function invcdf_gamma(a,x)
  y = ccall((:c_invcdf_gamma,lib), Float64, (Ref{Float64},Ref{Float64},), Float64(a), Float64(x))
  return y
end

function cdf_beta(a,b,x)
  y = ccall((:c_cdf_beta,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function pdf_beta(a,b,x)
  y = ccall((:c_pdf_beta,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function logpdf_beta(a,b,x)
  y = ccall((:c_logpdf_beta,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

function invcdf_beta(a,b,x)
  y = ccall((:c_invcdf_beta,lib), Float64, (Ref{Float64},Ref{Float64},Ref{Float64},), Float64(a), Float64(b), Float64(x))
  return y
end

end # end of module
