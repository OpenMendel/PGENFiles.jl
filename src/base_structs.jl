# TODO: DiffList

struct BitsVector{V} <: AbstractVector{UInt8}
    data::Base.RefValue{V}
    bits_per_element::UInt8
    size::UInt
end

function BitsVector(data::AbstractVector{UInt8}, bits_per_element, size)
    V = typeof(data)
    BitsVector{V}(Ref(data), bits_per_element, size)
end

struct ScatteredBitsVector{V} <: AbstractVector{UInt8}
    data_sectors::Vector{Base.RefValue{V}}
    bits_per_element::UInt8
    size::UInt
    n_blocks::UInt
end

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
    byte = x.data[][byte_index]
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

mutable struct DiffList{V,W,X,Y}
    len::UInt32
    sample_id_bases::Base.RefValue{V}
    last_component_sizes::Base.RefValue{W}
    has_genotypes::Bool
    genotypes::Union{BitsVector{X}, Nothing}
    sample_id_increments::Base.RefValue{Y}
end
