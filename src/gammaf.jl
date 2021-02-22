using SpecialFunctions
"""
    gammaf()
Variables: 0.5 <= x <= 5; 
Objective: Minimize gamma(x)
"""

function gammaf()
    best_known_objective = gamma(1.4616321)
    cplex_optimizer =
        JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)
    error_tolerances = [NaN64, 1e-3]
    # base partition holds the inflection points and the end-points
    milp_statistics = zeros(length(error_tolerances), 5)
    lp_statistics = zeros(length(error_tolerances), 5)

    base_partition = [0.5, 5] 
    # functions
    f = x -> gamma(x)

        # solve MILP for various error tolerances
        for i = 1:length(error_tolerances)
            error_tolerance = error_tolerances[i]
            # construct MILP relaxations
            milp = Model(cplex_optimizer)
            var, x_index, y_index = add_milp_relaxation(
                    f,
                    base_partition,
                    error_tolerance,
                    milp,
                    name = "1",
                )
            @objective(milp, Min, var[y_index])
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
            # construct LP relaxations
            lp = Model(cplex_optimizer)
            var, x_index, y_index = add_lp_relaxation(
                    f,
                    base_partition,
                    error_tolerance,
                    lp,
                    name = "1",
                )
            @objective(lp, Min, var[y_index])
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
    
        open("csv/gammaf_milp.csv", "w") do io
            writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
            writedlm(io, milp_statistics, ",")
        end
    
        open("csv/gammaf_lp.csv", "w") do io
            writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
            writedlm(io, lp_statistics, ",")
        end
    
end

