"""
    BitsVector{V}(data, bits_per_element, size)

Packed vector of 1, 2, or 4-bit entries. 
"""
struct BitsVector <: AbstractVector{UInt8}
    data::Ptr{UInt8}
    bits_per_element::UInt8
    size::UInt
end

"""
    ScatteredBitsVector{V}(data_sectors, bits_per_element, size, n_blocks)

Scattered packed vectors of 1, 2, or 4-bit entries. Length of each block is 2^16, 
except for the last one.
"""
struct ScatteredBitsVector{V} <: AbstractVector{UInt8}
    data_sectors::Vector{Base.RefValue{V}}
    bits_per_element::UInt8
    size::UInt
    n_blocks::UInt
end

"""
    ScatteredVector{T, V}(data_sectors, size, n_blocks)

Scattered vectors of type T. Length of each block is 2^16, except for the last one.
"""
struct ScatteredVector{T, V} <: AbstractVector{T}
    data_sectors::Vector{Base.RefValue{V}}
    size::Int
    n_blocks::Int
end

function Base.size(x::Union{BitsVector, ScatteredVector, ScatteredBitsVector})
    return (x.size,)
end

@inline function Base.getindex(x::BitsVector, i::Int)
    elements_per_byte = 8 ÷ x.bits_per_element
    byte_index = (i - 1) ÷ elements_per_byte + 1
    in_byte_index = (i - 1) % elements_per_byte
    byte = unsafe_load(x.data, byte_index)# x.data[][byte_index]
    (byte >> (x.bits_per_element * in_byte_index)) & mask_map[x.bits_per_element]
end

@inline function Base.getindex(x::ScatteredVector, i::Int)
    block_index = (i - 1) ÷ (2 ^ 16) + 1
    in_block_index = (i - 1) % (2 ^ 16) + 1
    x.data_sectors[block_index][][in_block_index]
end

@inline function Base.getindex(x::ScatteredBitsVector, i::Int)
    elements_per_byte = 8 ÷ x.bits_per_element
    block_index = (i - 1) ÷ (2 ^ 16) + 1
    in_block_index = ((i - 1) % (2 ^ 16)) ÷ elements_per_byte + 1
    in_byte_index = (i - 1) % elements_per_byte # 0-based
    byte = x.data_sectors[block_index][][in_block_index]
    (byte >> (x.bits_per_element * in_byte_index)) & mask_map[x.bits_per_element]
end

"""
    DiffList{V,W,X,Y}(len, sample_id_bases, last_component_sizes, has_genotypes, 
    genotypes, sample_id_increments)

Data structure for difflists.
"""
mutable struct DiffList{SIBT,GT<:Union{Nothing,BitsVector}}
    len::UInt32
    sample_id_bases::SIBT
    last_component_sizes::Ptr{UInt8}
    has_genotypes::Bool
    genotypes::GT
    sample_id_increments::Ptr{UInt8}
end