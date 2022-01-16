abstract type VariantIterator end

function offset_first_variant(x::Pgen)::UInt64
    x.header.variant_block_offsets[1]
end

@inline function Base.eltype(vi::VariantIterator)
    Variant
end

"""
    VariantIteratorBase(p::Pgen)
Variant iterator that iterates from the beginning of the Pgen file
"""
struct VariantIteratorBase <: VariantIterator
    p::Pgen
    v::Variant
end

function Base.iterate(vi::VariantIteratorBase,
    state=(1, offset_first_variant(vi.p)))
    if state[1] > vi.p.header.n_variants
        return nothing
    else
        idx = state[1]
        vi.v.index = state[1]
        vi.v.offset = state[2]
        vi.v.record_type = vi.p.header.variant_types[idx]
        vi.v.length = vi.p.header.variant_lengths[idx]
        nextstate = (idx + 1, state[2] + vi.v.length)
        return (vi.v, nextstate)
    end
end

@inline function Base.length(vi::VariantIteratorBase)
    vi.p.header.n_variants
end

@inline function Base.size(vi::VariantIteratorBase)
    (vi.p.header.n_variants, )
end

"""
    iterator(p::Pgen; startidx=1)
    
Retrieve a variant iterator for `p`.
"""
function iterator(p::Pgen; startidx=1)
    if startidx == 1
        v = Variant(0, 0, 0, 0)
        VariantIteratorBase(p, v)
    else
        @assert false "Not implemented."
    end
end

@inline function set_first_variant!(v::Variant, p::Pgen)
    v.index = 1
    v.offset = offset_first_variant(p)
    v.record_type = p.header.variant_types[1]
    v.length = p.header.variant_lengths[1]
    v
end
