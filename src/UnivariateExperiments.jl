module UnivariateExperiments

import Plots
import PolyhedralRelaxations

const PR = PolyhedralRelaxations

function check_package_api()
    formulation_data, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @info formulation_data
end

# This line changes the Plots backend to GR, which speeds up plotting.
# An alternative is `Plots.plotly()`, which uses Python's plotly library as the backend.
Plots.gr()

function plot(x::Vector{<:Real},y::Vector{<:Real},name::String)
    p = Plots.plot(x,y)
    Plots.savefig(p, name)
end

function plot_mip()
    xvals = collect(-1.0:0.01:1.0)
    yvals = [x^3 for x in xvals]
    plot(xvals,yvals,"plots/xpower3.pdf")
end

end # module
