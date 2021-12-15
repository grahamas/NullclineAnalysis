module NullclineAnalysis

using Simulation73
using OffsetArrays, StaticArrays
using StatsBase
using DataStructures: MutableLinkedList, ListNode, length
using NamedDims
using Contour: contour, lines
using IterTools
using ForwardDiff
using LinearAlgebra
using Base.Threads, LoopVectorization
using ThreadsX

include("types.jl")
export AbstractNullclineParams

include("stubs.jl")
export field_functions, phase_space_bounds, 
    derive_vector_fn, derive_vector_fn!,
    derive_jacobian_fn, derive_jacobian_fn!

include("fixedpoints.jl")
export calculate_fixedpoints, calculate_fixedpoints!

include("fixedpoint_stability.jl")
export fixedpoint_is_stable, fixedpoint_stability, 
    fixedpoint_is_oscillatory,
    count_stable_fps, filter_stable_fps,
    count_oscillatory_fps

end # module
