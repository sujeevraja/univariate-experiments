#!/usr/bin/env bash
#=
export JULIA_DEBUG=UnivariateExperiments
exec $HOME/bin/julia/julia-1.3.0/bin/julia --project="." --color=yes --startup-file=no -e \
    'include(popfirst!(ARGS))' "${BASH_SOURCE[0]}" "$@"
=#

using ArgParse
using UnivariateExperiments

s = ArgParseSettings()
@add_arg_table! s begin
    "--opt1"
        help = "an option with an argument"
    "--opt2", "-o"
        help = "another option with an argument"
        arg_type = Int
    "--flag1"
        help = "an option without an argument, i.e. a flag"
        action = :store_true
end

parsed_args = parse_args(ARGS, s)

@info parsed_args

UnivariateExperiments.greet()