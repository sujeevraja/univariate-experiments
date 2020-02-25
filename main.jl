#!/usr/bin/env bash
#=
export JULIA_DEBUG=UnivariateExperiments
exec $HOME/bin/julia/julia-1.3.0/bin/julia --project="." --color=yes --startup-file=no -e \
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

include("src/UnivariateExperiments.jl")
include("src/minlplib_helper.jl")
include("src/trig.jl")
include("src/ex4_1_1.jl")
include("src/fo7.jl")

s = ArgParseSettings()
@add_arg_table! s begin
    "--plot"
    help = "flag to generate the plots"
    action = :store_true
    "--minlplib"
    help = "flag to generate the results for minlplib"
    action = :store_true
    "--basic"
    help = "flag to collect the basic results for atomic functions"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)

if parsed_args["plot"]
    generate_plots()
end 

if parsed_args["minlplib"]
    @info "running the relaxations on trig"
    trig()
    @info "running the relaxations on ex4_1_1"
    ex4_1_1()
    @info "running the relaxations on fo7"
    fo7()
end 
# @info parsed_args

# UnivariateExperiments.check_package_api()
# plot()
