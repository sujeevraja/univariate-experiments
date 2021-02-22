#!/usr/bin/env bash
#=
export JULIA_DEBUG=UnivariateExperiments
exec julia --project="." --color=yes --startup-file=no -e \
    'include(popfirst!(ARGS))' "${BASH_SOURCE[0]}" "$@"
=#

using JuMP
using CPLEX
using PolyhedralRelaxations
using Printf
using Plots
using SparseArrays
using ArgParse
using DelimitedFiles

include("src/basics.jl")
include("src/plotting.jl")
include("src/minlplib_helper.jl")
include("src/trig.jl")
include("src/ex4_1_1.jl")
include("src/fo7.jl")
include("src/errorf.jl")
include("src/ramsey.jl")
include("src/gammaf.jl")

s = ArgParseSettings()
@add_arg_table! s begin
    "--plot", "-p"
    help = "generate plots"
    action = :store_true
    "--minlplib", "-m"
    help = "generate results for MINLPLIB instances"
    action = :store_false
    "--basic", "-b"
    help = "collect basic results for atomic functions"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

if parsed_args["plot"]
    generate_plots()
end

if parsed_args["minlplib"]
    # @info "running the relaxations on trig"
    # trig()
    # @info "running the relaxations on ex4_1_1"
    # ex4_1_1()
    # @info "running the relaxations on fo7"
    # fo7()
    @info "running the relaxations on errorf"
    errorf()
    @info "running the relaxations on ramsey"
    ramsey()
    @info "running the relaxations on gamma"
    gammaf()
end

if parsed_args["basic"]
    @info "generating error result table..."
    generate_error_table()
    @info "generated error resul table."
end
