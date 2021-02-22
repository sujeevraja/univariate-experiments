using SpecialFunctions
"""
    errorf()
``errorf'' is an NLP instance taken from 
https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir85-3206.pdf
The best known primal feasible objective value for this problem approximately 0. 
This example contains 5 nonlinear terms.
Variables: 0 <= x[1:2] <= 10; 
Objective: Minimize errorf(x[1] + x[2]) + sin(x[1]) * exp(-x[2]/10)
Constraints: x[1]^2 + x[2]^2 >= 10

The reformulation using multiple additional variables is as follows:
Variables: 
    0 <= x[1:2] <= 10
    x[3:4] - transformation variables
    -0.632121 <= x[3] <= 1.90484
    -1.90484 <= x[4] <= 0.632121
    y[1:8]
    0 <= y[1] <= 20
Constraints: 
    y[1] = x[1] + x[2]
    y[2] = errorf(y[1])
    x[3] = y[3] + y[4]
    x[4] = y[3] - y[4]
    y[3] = sin(x[1]) 
    y[4] = exp(-x[2]/10)
    y[5] = x[1]^2
    y[6] = x[2]^2 
    y[7] = x[3]^2
    y[8] = x[4]^2
Objective: Minimize y[2] + 0.25 * (y[7] - y[8])
"""

function errorf()
    best_known_objective = round(erf(3*π/2) + sin(3*π/2) * 1.0; digits=4)
    cplex_optimizer =
        JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)
    error_tolerances = [NaN64, 1e-3]
    # base partition holds the inflection points and the end-points
    milp_statistics = zeros(length(error_tolerances), 5)
    lp_statistics = zeros(length(error_tolerances), 5)

    base_partition = Dict{Int,Vector{Float64}}()
    base_partition[1] = [0.0, 1.0] # erf()
    base_partition[2] = [0.0, π, 2*π, 3*π, 10] # sin() 
    base_partition[3] = [0, 10.0] # exp
    base_partition[4] = [0, 10.0] # x[1]^2
    base_partition[5] = [0, 10.0] # x[2]^2
    base_partition[6] = [-0.632121, 1.90484] # x[3]^2 
    base_partition[7] = [-1.90484, 0.632121] # x[4]^2
    # functions
    functions = Dict{Int,Function}()
    functions[1] = x -> erf(x)
    functions[2] = x -> sin(x)
    functions[3] = x -> exp(-x/10)
    functions[4] = x -> x^2
    functions[5] = x -> x^2
    functions[6] = x -> x^2
    functions[7] = x -> x^2

        # solve MILP for various error tolerances
        for i = 1:length(error_tolerances)
            error_tolerance = error_tolerances[i]
            milp_var = Dict{Int,Any}()
            x_indexes = Dict{Int,Int}()
            y_indexes = Dict{Int,Int}()
            # construct MILP relaxations
            milp = Model(cplex_optimizer)
            x_lb = [0, 0, -0.632121, -1.90484]
            x_ub = [10, 10, 1.90484, 0.632121]
            @variable(milp, x_lb[i] <= x[i=1:4] <= x_ub[i])
            @variable(milp, y[1:8])
            set_lower_bound(y[1], 0)
            set_upper_bound(y[1], 20)
            @constraint(milp, y[1] == x[1] + x[2])
            @constraint(milp, x[3] == y[3] + y[4])
            @constraint(milp, x[4] == y[3] - y[4])
            @objective(milp, Min, y[2] + 0.25 * (y[7] - y[8]))
            var_map = Dict{Int,Any}()
            var_map[1] = (y[1], y[2]) # (indep, dep)
            var_map[2] = (x[1], y[3])
            var_map[3] = (x[2], y[4])
            var_map[4] = (x[1], y[5])
            var_map[5] = (x[2], y[6])
            var_map[6] = (x[3], y[7])
            var_map[7] = (x[4], y[8])
            for j = 1:7
                var, x_index, y_index = add_milp_relaxation(
                    functions[j],
                    base_partition[j],
                    error_tolerance,
                    milp,
                    name = string(j),
                )
                milp_var[j] = var
                x_indexes[j] = x_index
                y_indexes[j] = y_index
                @constraint(milp, var_map[j][1] == var[x_index])
                @constraint(milp, var_map[j][2] == var[y_index])
            end
            _, solve_time, solve_bytes_alloc, sec_in_gc = @timed optimize!(milp)
            relaxation_objective = objective_value(milp)
            @assert relaxation_objective ≤ best_known_objective
            relative_gap =
                abs(best_known_objective - relaxation_objective) / (abs(relaxation_objective) + 1e-4)
            milp_statistics[i, :] = [
                error_tolerance,
                best_known_objective,
                relaxation_objective,
                relative_gap,
                round(solve_time, digits = 2),
            ]
        end

        for i = 1:length(error_tolerances)
            error_tolerance = error_tolerances[i]
            lp_var = Dict{Int,Any}()
            x_indexes = Dict{Int,Int}()
            y_indexes = Dict{Int,Int}()
            # construct LP relaxations
            lp = Model(cplex_optimizer)
            x_lb = [0, 0, -0.632121, -1.90484]
            x_ub = [10, 10, 1.90484, 0.632121]
            @variable(lp, x_lb[i] <= x[i=1:4] <= x_ub[i])
            @variable(lp, y[1:8])
            set_lower_bound(y[1], 0)
            set_upper_bound(y[1], 20)
            @constraint(lp, y[1] == x[1] + x[2])
            @constraint(lp, x[3] == y[3] + y[4])
            @constraint(lp, x[4] == y[3] - y[4])
            @objective(lp, Min, y[2] + 0.25 * (y[7] - y[8]))
            var_map = Dict{Int,Any}()
            var_map[1] = (y[1], y[2]) # (indep, dep)
            var_map[2] = (x[1], y[3])
            var_map[3] = (x[2], y[4])
            var_map[4] = (x[1], y[5])
            var_map[5] = (x[2], y[6])
            var_map[6] = (x[3], y[7])
            var_map[7] = (x[4], y[8])
            for j = 1:7
                var, x_index, y_index = add_lp_relaxation(
                    functions[j],
                    base_partition[j],
                    error_tolerance,
                    lp,
                    name = string(j),
                )
                lp_var[j] = var
                x_indexes[j] = x_index
                y_indexes[j] = y_index
                @constraint(lp, var_map[j][1] == var[x_index])
                @constraint(lp, var_map[j][2] == var[y_index])
            end
            _, solve_time, solve_bytes_alloc, sec_in_gc = @timed optimize!(lp)
            relaxation_objective = objective_value(lp)
            @assert relaxation_objective ≤ best_known_objective
            relative_gap =
                abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
            lp_statistics[i, :] = [
                error_tolerance,
                best_known_objective,
                relaxation_objective,
                relative_gap,
                round(solve_time, digits = 2),
            ]
        end
    
        open("csv/erf_milp.csv", "w") do io
            writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
            writedlm(io, milp_statistics, ",")
        end
    
        open("csv/erf_lp.csv", "w") do io
            writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
            writedlm(io, lp_statistics, ",")
        end
    
end

