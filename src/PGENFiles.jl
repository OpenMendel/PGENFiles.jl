module PGENFiles
using BitIntegers
import Mmap: mmap
import Base: unsafe_load
export Pgen, iterator, n_samples, n_variants, get_genotypes, get_genotypes!
export alt_allele_dosage, alt_allele_dosage!, ref_allele_dosage, ref_allele_dosage!
export write_PGEN
BitIntegers.@define_integers 24
const variant_type_lengths = Dict(
    0x00 => (4, 1), 0x01 => (4, 2), 0x02 => (4, 3), 0x03 => (4, 4),
    0x04 => (8, 1), 0x05 => (8, 2), 0x06 => (8, 3), 0x07 => (8, 4)
)
const bytes_to_UInt = Dict(0x01 => UInt8, 0x02 => UInt16, 0x03 => UInt24, 0x04 => UInt32, 0x08 => UInt64)
const mask_map = [0x01, 0x03, 0x00, 0x0f, 0x00, 0x00, 0x00, 0xff]

@inline function Base.unsafe_load(p::Ptr{UInt24}, i::Int=1)
    p_UInt8 = reinterpret(Ptr{UInt8}, p)
    p_target = p_UInt8 + (i - 1) * 3
    a = unsafe_wrap(Array{UInt8}, p_target, (3,))
    reinterpret(UInt24, a)[1]
end

@inline ceil_int(x::Integer, y::Integer) = (x ÷ y) + (x % y != 0)
include("uleb128.jl")
include("internal_structs.jl")
include("structs.jl")
include("header.jl")
include("difflist.jl")
include("iterator.jl")
include("genotype.jl")
include("dosage.jl")
include("write.jl")
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)
end
