module UnivariateExperiments

using PolyhedralRelaxations

function greet()
    formulation_data, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @info formulation_data
end

end # module
