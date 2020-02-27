import PolyhedralRelaxations

const PR = PolyhedralRelaxations

struct InputData
    f::Function
    f_name::String
    base_partition::Vector{Float64}
end

struct ErrorResults
    num_vars_for_errs::Vector{Int64}
    errs_for_num_vars::Vector{Float64}
end

function get_error_results(input_data::InputData,err_tols::Vector{Float64},num_vars::Vector{Int64})::ErrorResults
    num_vars_for_errs = Int64[]
    err_var_dict = Dict{Float64,Int64}
    for err_tol in err_tols
        milp_data, _ = PR.construct_milp_relaxation(input_data.f,
            input_data.base_partition, error_tolerance=err_tol)
        push!(num_vars_for_errs, PR.get_num_binary_variables(milp_data))
    end

    errs_for_num_vars = Float64[]
    for num in num_vars
        milp_data, _ = PR.construct_milp_relaxation(input_data.f, input_data.base_partition,
            num_additional_binary_variables=num)
        push!(errs_for_num_vars, PR.get_error_bound(milp_data))
    end

    return ErrorResults(num_vars_for_errs, errs_for_num_vars)
end

function generate_error_table()
    inputs = [
        InputData(sin, "sin", [0,π,2*π]),
        InputData(x -> x^3, "x^3", [-1.0,0.0,1.0]),
        InputData(x -> x * abs(x), "x|x|", [-2.0,0.0,2.0]),
        InputData(x -> 1 / (1+exp(-x)), "logistic(x)", [-5.0,0.0,5.0])
    ]
    err_tols = [NaN64,0.1,0.01]
    num_vars = [0,50,100]
    results = ErrorResults[]
    for input_data in inputs
        push!(results, get_error_results(input_data,err_tols,num_vars))
    end

    headers = [
        "name",
        "base_partition",
        "num_partitions_error_0",
        "num_partitions_error_0.1",
        "num_partitions_error_0.01",
        "strength_0_vars",
        "strength_50_vars",
        "strength_100_vars",
    ]

    open("out/error_results.csv", "w") do io
        write(io, get_string_from_list(headers))
        write(io, "\n")
        for i in 1:length(results)
            input_data = inputs[i]
            error_results = results[i]

            row = String[]
            push!(row, input_data.f_name)
            push!(row, string(input_data.base_partition))
            for j in 1:length(err_tols)
                push!(row, string(error_results.num_vars_for_errs[j]))
            end

            for j in 1:length(num_vars)
                push!(row, string(error_results.errs_for_num_vars[j]))
            end

            write(io, get_string_from_list(row))
            write(io, "\n")
        end
    end
end

function get_string_from_list(l::Vector{String})::String
    s = l[1]
    for i in 2:length(l)
        s = s * ";"
        s = s * l[i]
    end
    return s
end
