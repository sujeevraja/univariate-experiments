import Cbc
import Plots
import PolyhedralRelaxations
using JuMP
using SparseArrays

const PR = PolyhedralRelaxations

# This line changes the Plots backend to GR, which speeds up plotting.
# An alternative is `Plots.plotly()`, which uses Python's plotly library as the backend.
Plots.gr()

base_partition = collect(-1.0:0.25:1.0)

function get_mip_values()::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    cbc_opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    milp = Model(cbc_opt)
    f = x->x^3
    milp_data, function_data = PR.construct_milp_relaxation(f,base_partition)
    lb, ub = PR.get_variable_bounds(milp_data)

    @info "x bounds: $(lb[milp_data.x_index]), $(ub[milp_data.x_index])"
    @info "y bounds: $(lb[milp_data.y_index]), $(ub[milp_data.y_index])"

    # Add variables
    num_vars = PR.get_num_variables(milp_data)
    binary_indices = PR.get_binary_indices(milp_data)
    @variable(milp, lb[i] <= x[i=1:num_vars] <= ub[i])
    for k in binary_indices
        set_binary(x[k])
    end
    for k in 1:num_vars
        set_name(x[k], PR.get_variable_names(milp_data)[k])
    end

    # Add constraints
    A, b = PR.get_eq_constraint_matrices(milp_data)
    @constraint(milp, A * x .== b)
    A, b = PR.get_leq_constraint_matrices(milp_data)
    @constraint(milp, A * x .<= b)

    # Add objective.
    set_objective_function(milp, x[milp_data.y_index])

    # Collect solutions from different linear objectives.
    xs = collect(-1.0:0.01:1.0)
    max_ys = Float64[]
    min_ys = Float64[]

    for p in xs
        set_lower_bound(x[milp_data.x_index], p)
        set_upper_bound(x[milp_data.x_index], p)

        set_objective_sense(milp, MOI.MAX_SENSE)
        optimize!(milp)
        @assert termination_status(milp) == MOI.OPTIMAL
        yval = value.(x[milp_data.y_index])
        obj = objective_value(milp)
        @assert isapprox(yval, obj, atol=1e-5) "obj $obj and y $yval different in max problem"
        push!(max_ys, yval)

        set_objective_sense(milp, MOI.MIN_SENSE)
        optimize!(milp)
        @assert termination_status(milp) == MOI.OPTIMAL
        yval = value.(x[milp_data.y_index])
        obj = objective_value(milp)
        @assert isapprox(yval, obj, atol=1e-5) "obj $obj and y $yval different in min problem"
        push!(min_ys, yval)
    end

    return xs,max_ys,min_ys
end

function plot_mip()
    xs,max_ys,min_ys = get_mip_values()
    fs = [x^3 for x in xs]
    p = Plots.plot(xs,hcat(max_ys,fs,min_ys), legend=false)
    Plots.savefig(p, "out/mip.pdf")
end

function plot_lp()
    f = x -> x^3
    xs = collect(-1.0:0.01:1.0)
    fs = [x^3 for x in xs]
    p = Plots.plot(xs,fs,legend=false,color=:orange)

    lp_data, fn_data = PR.construct_lp_relaxation(f, base_partition)
    cbc_opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    lp = Model(cbc_opt)
    lb, ub = PR.get_variable_bounds(lp_data)
    @info "x bounds: $(lb[lp_data.x_index]), $(ub[lp_data.x_index])"
    @info "y bounds: $(lb[lp_data.y_index]), $(ub[lp_data.y_index])"

    num_vars = PR.get_num_variables(lp_data)
    @variable(lp, lb[i] <= y[i=1:num_vars] <= ub[i])
    A, b = PR.get_eq_constraint_matrices(lp_data)
    @constraint(lp, A * y .== b)

    for k in 1:num_vars
        set_name(y[k], PR.get_variable_names(lp_data)[k])
    end

    y_min, y_max = -1.0,1.0

    res = π/100
    α_min = res
    α_max = (π/2) - res
    for α ∈ collect(α_min:res:α_max)
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
        Plots.plot!(line_xs, line_ys, color=:blue)

        @objective(lp, Min, yvar - (xvar * tan(α)))
        optimize!(lp)
        @assert termination_status(lp) == MOI.OPTIMAL
        c = objective_value(lp)
        x_min = max((y_min - c) / m, -1.0)
        x_max = min((y_max - c) / m, 1.0)
        line_xs = [x_min, x_max]
        line_ys = [(m*x) + c for x in line_xs]
        Plots.plot!(line_xs, line_ys, color=:green)
    end

    Plots.savefig(p, "out/lp.pdf")
end
