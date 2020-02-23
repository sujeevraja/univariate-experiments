import Cbc
import Plots
import PolyhedralRelaxations
using Printf
using JuMP
using SparseArrays

const PR = PolyhedralRelaxations

# This line changes the Plots backend to GR, which speeds up plotting.
# An alternative is `Plots.plotly()`, which uses Python's plotly library as the backend.
Plots.gr()

function get_mip_values()::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    f = x->x^3
    f_dash = x->3*(x^2)
    part = collect(-1.0:0.5:1.0)
    formulation_data, function_data = PR.construct_milp_relaxation(f,f_dash,part)
    cbc_opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    milp = Model(cbc_opt)
    lb, ub = PR.get_variable_bounds(formulation_data)

    @info "x bounds: $(lb[formulation_data.x_index]), $(ub[formulation_data.x_index])"
    @info "y bounds: $(lb[formulation_data.y_index]), $(ub[formulation_data.y_index])"

    # Add variables
    num_vars = PR.get_num_variables(formulation_data)
    binary_indices, _ = findnz(formulation_data.binary)
    @variable(milp, lb[i] <= x[i=1:num_vars] <= ub[i])
    for k in binary_indices
        set_binary(x[k])
    end
    for k in 1:num_vars
        set_name(x[k], PR.get_variable_names(formulation_data)[k])
    end

    # Add constraints
    A, b = PR.get_eq_constraint_matrices(formulation_data)
    @constraint(milp, A * x .== b)
    A, b = PR.get_leq_constraint_matrices(formulation_data)
    @constraint(milp, A * x .<= b)

    # Add objective.
    set_objective_function(milp, x[formulation_data.y_index])

    # Collect solutions from different linear objectives.
    xs = collect(-1.0:0.01:1.0)
    max_ys = Float64[]
    min_ys = Float64[]

    for p in xs
        set_lower_bound(x[formulation_data.x_index], p)
        set_upper_bound(x[formulation_data.x_index], p)

        set_objective_sense(milp, MOI.MAX_SENSE)
        optimize!(milp)
        @assert termination_status(milp) == MOI.OPTIMAL
        yval = value.(x[formulation_data.y_index])
        obj = objective_value(milp)
        @assert isapprox(yval, obj, atol=1e-5) "obj $obj and y $yval different in max problem"
        push!(max_ys, yval)
        # @printf "x %.3f y %.3f max obj %.3f\n" p yval objective_value(milp)

        set_objective_sense(milp, MOI.MIN_SENSE)
        optimize!(milp)
        @assert termination_status(milp) == MOI.OPTIMAL
        yval = value.(x[formulation_data.y_index])
        obj = objective_value(milp)
        @assert isapprox(yval, obj, atol=1e-5) "obj $obj and y $yval different in min problem"
        push!(min_ys, yval)
    end

    return xs,max_ys,min_ys
end

function plot_mip()
    xs,max_ys,min_ys = get_mip_values()
    fs = [x^3 for x in xs]
    p = Plots.plot(xs,hcat(max_ys,fs,min_ys))
    Plots.savefig(p, "out/mip.pdf")
end

function plot_lp()
    f = x -> x^3
    xs = collect(-1.0:0.01:1.0)
    fs = [x^3 for x in xs]
    p = Plots.plot(xs,fs,legend=false)

    cbc_opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    lp_data, fn_data = PR.construct_lp_relaxation(f, collect(-1.0:0.25:1.0))
    lp = Model(cbc_opt)
    lb, ub = PR.get_variable_bounds(lp_data)
    num_vars = PR.get_num_variables(lp_data)
    @variable(lp, lb[i] <= y[i=1:num_vars] <= ub[i])
    A, b = PR.get_eq_constraint_matrices(lp_data)
    @constraint(lp, A * y .== b)

    y_min, y_max = -1.0,1.0

    α_max = π/2
    res = π/100
    for α ∈ collect(res:res:(α_max-res))
        m = tan(α)
        xvar = y[lp_data.x_index]
        yvar = y[lp_data.y_index]

        @objective(lp, Max, yvar - (xvar * tan(α)))
        optimize!(lp)
        @assert termination_status(lp) == MOI.OPTIMAL
        c = objective_value(lp)
        x_min = max((y_min - c) / m, -1.0)
        x_max = min((y_max - c) / m, 1.0)
        line_xs = [x_min, x_max]
        line_ys = [(m*x) + c for x in line_xs]
        Plots.plot!(line_xs, line_ys)

        @objective(lp, Min, yvar - (xvar * tan(α)))
        optimize!(lp)
        @assert termination_status(lp) == MOI.OPTIMAL
        c = objective_value(lp)
        x_min = max((y_min - c) / m, -1.0)
        x_max = min((y_max - c) / m, 1.0)
        line_xs = [x_min, x_max]
        line_ys = [(m*x) + c for x in line_xs]
        Plots.plot!(line_xs, line_ys)
    end

    Plots.savefig(p, "out/lp.pdf")
end
