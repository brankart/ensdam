push!(LOAD_PATH, "./")
import ensdam
y=ensdam.cdf_gaussian(1.)
println(y)
y=ensdam.cdf_gamma(0.5,0.5)
println((y+1)*0.5)
