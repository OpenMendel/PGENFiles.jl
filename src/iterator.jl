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
end

function Base.iterate(vi::VariantIteratorBase,
    state=(1, offset_first_variant(vi.p)))
    if state[1] > vi.p.header.n_variants
        return nothing
    else
        idx = state[1]
        v = Variant(idx, state[2], vi.p.header.variant_types[idx], vi.p.header.variant_lengths[idx])
        nextstate = (idx + 1, state[2] + vi.p.header.variant_lengths[idx])
        return (v, nextstate)
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
        VariantIteratorBase(p)
    else
        @error "Not implemented."
    end
end