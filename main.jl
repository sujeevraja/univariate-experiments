#!/usr/bin/env bash
#=
JULIA_DEBUG=UnivariateExperiments
export JULIA_DEBUG
exec $HOME/bin/julia/julia-1.3.0/bin/julia --project="." --color=yes --startup-file=no -e \
    'include(popfirst!(ARGS))' "${BASH_SOURCE[0]}" "$@"
=#

@show ARGS

using UnivariateExperiments

UnivariateExperiments.greet()
