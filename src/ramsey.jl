"""
    ramsey()
``ramsey'' is an MINLPLib test instance with binary and continuous variables
(see http://www.minlplib.org/ramsey.html). The best known primal feasible
objective value for this problem is 2.48746864. This example contains 22 nonlinear
terms, 11 are x**0.25 and remaining are log(x)
"""

function post_ramsey_linear_constraints!(m)
@variable(m, objvar)
    x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
    @variable(m, x[x_Idx])
    setlowerbound(x[1], 3.0)
    setupperbound(x[1], 3.0)
    setlowerbound(x[2], 3.0)
    setlowerbound(x[3], 3.0)
    setlowerbound(x[4], 3.0)
    setlowerbound(x[5], 3.0)
    setlowerbound(x[6], 3.0)
    setlowerbound(x[7], 3.0)
    setlowerbound(x[8], 3.0)
    setlowerbound(x[9], 3.0)
    setlowerbound(x[10], 3.0)
    setlowerbound(x[11], 3.0)
    setlowerbound(x[12], 0.95)
    setlowerbound(x[13], 0.95)
    setlowerbound(x[14], 0.95)
    setlowerbound(x[15], 0.95)
    setlowerbound(x[16], 0.95)
    setlowerbound(x[17], 0.95)
    setlowerbound(x[18], 0.95)
    setlowerbound(x[19], 0.95)
    setlowerbound(x[20], 0.95)
    setlowerbound(x[21], 0.95)
    setlowerbound(x[22], 0.95)
    setlowerbound(x[23], 0.05)
    setupperbound(x[23], 0.05)
    setlowerbound(x[24], 0.05)
    setupperbound(x[24], 0.0575)
    setlowerbound(x[25], 0.05)
    setupperbound(x[25], 0.066125)
    setlowerbound(x[26], 0.05)
    setupperbound(x[26], 0.07604375)
    setlowerbound(x[27], 0.05)
    setupperbound(x[27], 0.0874503125)
    setlowerbound(x[28], 0.05)
    setupperbound(x[28], 0.100567859375)
    setlowerbound(x[29], 0.05)
    setupperbound(x[29], 0.11565303828125)
    setlowerbound(x[30], 0.05)
    setupperbound(x[30], 0.133000994023437)
    setlowerbound(x[31], 0.05)
    setupperbound(x[31], 0.152951143126953)
    setlowerbound(x[32], 0.05)
    setupperbound(x[32], 0.175893814595996)
    setlowerbound(x[33], 0.05)
    setupperbound(x[33], 0.202277886785395)


    # ----- Constraints ----- #
    @constraint(m, e12, -x[1]+x[2]-x[23] == 0.0)
    @constraint(m, e13, -x[2]+x[3]-x[24] == 0.0)
    @constraint(m, e14, -x[3]+x[4]-x[25] == 0.0)
    @constraint(m, e15, -x[4]+x[5]-x[26] == 0.0)
    @constraint(m, e16, -x[5]+x[6]-x[27] == 0.0)
    @constraint(m, e17, -x[6]+x[7]-x[28] == 0.0)
    @constraint(m, e18, -x[7]+x[8]-x[29] == 0.0)
    @constraint(m, e19, -x[8]+x[9]-x[30] == 0.0)
    @constraint(m, e20, -x[9]+x[10]-x[31] == 0.0)
    @constraint(m, e21, -x[10]+x[11]-x[32] == 0.0)
    @constraint(m, e22, 0.03*x[11]-x[33] <= 0.0)

    println("creating lifted variables and constraints")
    @variable(m, y[1:22])
    @constraint(m, y[1] == 3^0.25)
    @constraint(m, e1, 0.759835685651593*y[1]-x[12]-x[23] == 0.0)
    @constraint(m, e2, 0.77686866556676*y[2]-x[13]-x[24] == 0.0)
    @constraint(m, e3, 0.794283468039448*y[3]-x[14]-x[25] == 0.0)
    @constraint(m, e4, 0.812088652256959*y[4]-x[15]-x[26] == 0.0)
    @constraint(m, e5, 0.830292969275008*y[5]-x[16]-x[27] == 0.0)
    @constraint(m, e6, 0.848905366318769*y[6]-x[17]-x[28] == 0.0)
    @constraint(m, e7, 0.867934991180342*y[7]-x[18]-x[29] == 0.0)
    @constraint(m, e8, 0.88739119671479*y[8]-x[19]-x[30] == 0.0)
    @constraint(m, e9, 0.907283545436972*y[9]-x[20]-x[31] == 0.0)
    @constraint(m, e10, 0.92762181422141*y[10]-x[21]-x[32] == 0.0)
    @constraint(m, e11, 0.948415999107521*y[11]-x[22]-x[33] == 0.0)
    @constraint(m, e23, -(0.95*y[12]+0.9025*y[13]+0.857375*y[14]+0.81450625*y[15]+0.7737809375*y[16]+0.735091890625*y[17]+0.69833729609375*y[18]+0.663420431289062*y[19]+0.630249409724609*y[20]+0.598736939238379*y[21]+11.3760018455292*y[22])-objvar == 0.0)
    
    # ----- Objective ----- #
    @objective(m, Min, objvar)

    y_indexes_1 = collect(2:11)
    y_indexes_2 = collect(12:22)
    x_indexes_1 = collect(2:11)
    x_indexes_2 = collect(12:22)
    partitions_1 = []
    partitions_2 = []
    for index in x_indexes_1 
        ub = 1e2
        (has_upper_bound(x[index])) && (ub = upper_bound(x[index]))
        push!(partitions_1, [lower_bound(x[index]), ub])
    end 
    for index in x_indexes_2 
        ub = 1e2
        (has_upper_bound(x[index])) && (ub = upper_bound(x[index]))
        push!(partitions_2, [lower_bound(x[index]), ub])
    end 

    return y_indexes_1, x_indexes_1, partitions_1, y_indexes_2, x_indexes_2, partitions_2
end 

function ramsey()
    best_known_objective = -2.48746864
    cplex_optimizer =
        JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0)
    error_tolerances = [NaN64, 1e-1, 1e-2]
    # base partition holds the inflection points and the end-points for the trigonometric functions
    milp_statistics = zeros(length(error_tolerances), 5)
    lp_statistics = zeros(length(error_tolerances), 5)
    func_1 = x -> x^0.25
    func_2 = x -> log(x)
    # solve MILP for various error tolerances
    for i = 1:length(error_tolerances)
        error_tolerance = error_tolerances[i]
        milp = Model(cplex_optimizer)
        y_indexes_1, x_indexes_1, base_partition_1, 
        y_indexes_2, x_indexes_2, base_partition_2 = post_ramsey_linear_constraints!(milp)
        # construct MILP relaxations
        for j = 1:length(y_indexes_1)
            var, x_index, y_index = add_milp_relaxation(
                func_1,
                base_partition_1[j],
                error_tolerance,
                milp,
                name = "1_" * string(j),
            )
            @constraint(milp, var[x_index] == milp[:x][x_indexes_1[j]])
            @constraint(milp, var[y_index] == milp[:y][y_indexes_1[j]])
        end
        for j = 1:length(y_indexes_2)
            var, x_index, y_index = add_milp_relaxation(
                func_2,
                base_partition_2[j],
                error_tolerance,
                milp,
                name = "2_" * string(j),
            )
            @constraint(milp, var[x_index] == milp[:x][x_indexes_2[j]])
            @constraint(milp, var[y_index] == milp[:y][y_indexes_2[j]])
        end
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
        y_indexes_1, x_indexes_1, base_partition_1, 
        y_indexes_2, x_indexes_2, base_partition_2 = post_ramsey_linear_constraints!(lp)
        # construct MILP relaxations
        for j = 1:length(y_indexes_1)
            var, x_index, y_index = add_lp_relaxation(
                func_1,
                base_partition_1[j],
                error_tolerance,
                lp,
                name = "1_" * string(j),
            )
            @constraint(lp, var[x_index] == lp[:x][x_indexes_1[j]])
            @constraint(lp, var[y_index] == lp[:y][y_indexes_1[j]])
        end
        for j = 1:length(y_indexes_2)
            var, x_index, y_index = add_lp_relaxation(
                func_2,
                base_partition_2[j],
                error_tolerance,
                lp,
                name = "2_" * string(j),
            )
            @constraint(lp, var[x_index] == lp[:x][x_indexes_2[j]])
            @constraint(lp, var[y_index] == lp[:y][y_indexes_2[j]])
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

    open("csv/ramsey_milp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, milp_statistics, ",")
    end

    open("csv/ramsey_lp.csv", "w") do io
        writedlm(io, ["error_tol" "g_opt" "relaxation" "relative_gap" "time"], ',')
        writedlm(io, lp_statistics, ",")
    end

end 