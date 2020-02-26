"""
    ex4_1_1(; formulation=1)
``trig'' is an MINLPLib test instance with one variable and power terms
(see http://www.minlplib.org/ex4_1_1.html). The global objective value for 
this problem is -7.48731237. This code illustrates the strength of the 
relaxation obtained by introducing one artificial variable for each non-linear
term as against one variable for the the entire univariate function.
Variables: -2 <= x <= 11; start = 10
Objective: Minimize 0.1 + x^6 - 2.08 * x^5 + 0.4875 * x^4 + 7.1 * x^3 - 3.95 * x^2 - x

The reformulation using one additional variable is as follows:
Variables:
    -2 <= x <= 11; start = 10
    y 
Constraints: y = 0.1 + x^6 - 2.08 * x^5 + 0.4875 * x^4 + 7.1 * x^3 - 3.95 * x^2 - x
Objective: Minimize y 

The reformulation using multiple additional variables is as follows:
Variables: 
    -2 <= x <= 11; start = 10
    y[1:5]
Constraints: 
    y[1] = x^6
    y[2] = -2.08 * x^5 
    y[3] = 0.4875 * x^4
    y[4] = 7.1 * x^3
    y[5] = -3.95 * x^2
Objective: Minimize 0.1 + y[1] + y[2] + y[3] + y[4] + y[5] - x 
"""

function ex4_1_1(; formulation = 1)
    if formulation == 1
        ex4_1_1_single()
    else
        ex4_1_1_multiple()
    end
end

function ex4_1_1_single()
    best_known_objective = -7.48731237
    cplex_optimizer =
        JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)
    error_tolerances = [NaN64, 1e-1, 1e-2]
    base_partition = [-2.0, -1.7747, 0.1753, 11.0]
    milp_statistics = zeros(length(error_tolerances), 5)
    lp_statistics = zeros(length(error_tolerances), 5)
    func = x -> 0.1 + x^6 - 2.08 * x^5 + 0.4875 * x^4 + 7.1 * x^3 - 3.95 * x^2 - x
    for i = 1:length(error_tolerances)
        error_tolerance = error_tolerances[i]
        # construct MILP relaxations
        milp = Model(cplex_optimizer)
        milp_var, x_index, y_index =
            add_milp_relaxation(func, base_partition, error_tolerance, milp)
        @objective(milp, Min, milp_var[y_index])
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed optimize!(milp)
        relaxation_objective = objective_value(milp)
        @assert relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
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
        lp_var, x_index, y_index =
            add_lp_relaxation(func, base_partition, error_tolerance, lp)
        @objective(lp, Min, lp_var[y_index])
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

    open("out/single_ex4_1_1_milp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, milp_statistics, ",")
    end

    open("out/single_ex4_1_1_lp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, lp_statistics, ",")
    end
end

function ex4_1_1_multiple()
    best_known_objective = -7.48731237
    cplex_optimizer =
        JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)
    error_tolerances = [NaN64, 1e-1, 1e-2]
    base_partition = [-2.0, 0.0, 11.0]
    milp_statistics = zeros(length(error_tolerances), 5)
    lp_statistics = zeros(length(error_tolerances), 5)
    functions = Dict{Int,Function}()
    functions[1] = x -> x^6
    functions[2] = x -> -2.08 * x^5
    functions[3] = x -> 0.4875 * x^4
    functions[4] = x -> 7.1 * x^3
    functions[5] = x -> -3.95 * x^2
    for i = 1:length(error_tolerances)
        error_tolerance = error_tolerances[i]
        milp = Model(cplex_optimizer)
        milp_var = Dict{Int,Any}()
        x_indexes = Dict{Int,Int}()
        y_indexes = Dict{Int,Int}()
        # construct MILP relaxations
        for j = 1:5
            var, x_index, y_index = add_milp_relaxation(
                functions[j],
                base_partition,
                error_tolerance,
                milp,
                name = string(j),
            )
            milp_var[j] = var
            x_indexes[j] = x_index
            y_indexes[j] = y_index
            (j == 1) && (continue)
            @constraint(milp, var[x_index] == milp_var[1][x_indexes[1]])
        end
        @objective(
            milp,
            Min,
            0.1 + sum(milp_var[j][y_indexes[j]] for j = 1:5) - milp_var[1][x_indexes[1]]
        )
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed optimize!(milp)
        relaxation_objective = objective_value(milp)
        @assert relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
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
        lp = Model(cplex_optimizer)
        lp_var = Dict{Int,Any}()
        x_indexes = Dict{Int,Int}()
        y_indexes = Dict{Int,Int}()
        # construct LP relaxations
        for j = 1:5
            var, x_index, y_index = add_lp_relaxation(
                functions[j],
                base_partition,
                error_tolerance,
                lp,
                name = string(j),
            )
            lp_var[j] = var
            x_indexes[j] = x_index
            y_indexes[j] = y_index
            (j == 1) && (continue)
            @constraint(lp, var[x_index] == lp_var[1][x_indexes[1]])
        end
        @objective(
            lp,
            Min,
            0.1 + sum(lp_var[j][y_indexes[j]] for j = 1:5) - lp_var[1][x_indexes[1]]
        )
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

    open("out/multiple_ex4_1_1_milp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, milp_statistics, ",")
    end

    open("out/multiple_ex4_1_1_lp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, lp_statistics, ",")
    end
end
