"""
    PgenVariantIterator(p::Pgen)
Variant iterator that iterates from the beginning of the Pgen file
"""
struct PgenVariantIterator <: VariantIterator
    p::Pgen
    v::Variant
end

function offset_first_variant(x::Pgen)::UInt64
    x.header.variant_block_offsets[1]
end

@inline function Base.eltype(vi::PgenVariantIterator)
    Variant
end

function Base.iterate(vi::PgenVariantIterator,
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

@inline function Base.length(vi::PgenVariantIterator)
    vi.p.header.n_variants
end

@inline function Base.size(vi::PgenVariantIterator)
    (vi.p.header.n_variants, )
end

# Overload collect for PgenVariantIterator
function Base.collect(vi::PgenVariantIterator)
    result = Vector{Variant}(undef, length(vi))
    i = 1
    for v in vi
        result[i] = PgenVariant(v.index, v.offset, v.record_type, v.length)
        i += 1
    end
    return result
end

"""
    iterator(p::Pgen; startidx=1)
    
Retrieve a variant iterator for `p`.
"""
function iterator(p::Pgen; startidx=1)
    if startidx == 1
        v = PgenVariant(0, 0, 0, 0)
        PgenVariantIterator(p, v)
    else
        @assert false "Not implemented."
    end
end

@inline function set_first_variant!(v::PgenVariant, p::Pgen)
    v.index = 1
    v.offset = offset_first_variant(p)
    v.record_type = p.header.variant_types[1]
    v.length = p.header.variant_lengths[1]
    v
end

@inline function chrom(p::Pgen, v::PgenVariant)
    return string(p.pvar_df[v.index, Symbol("#CHROM")])
end

@inline function pos(p::Pgen, v::PgenVariant)
    return p.pvar_df[v.index, Symbol("POS")]
end

@inline function rsid(p::Pgen, v::PgenVariant)
    return string(p.pvar_df[v.index, Symbol("ID")])
end

function alleles(p::Pgen, v::PgenVariant)
    return [p.pvar_df[v.index, Symbol("REF")], p.pvar_df[v.index, Symbol("ALT")]]
end

function alt_allele(p::Pgen, v::PgenVariant)
    return p.pvar_df[v.index, Symbol("ALT")]
end

function ref_allele(p::Pgen, v::PgenVariant)
    return p.pvar_df[v.index, Symbol("REF")]
end

function load_values!(arr::AbstractArray, p::Pgen, v::PgenVariant; genobuf = Vector{UInt8}(undef, n_samples), genoldbuf=nothing)
    alt_allele_dosage!(arr, genobuf, p, v; genoldbuf=genoldbuf)
    arr
end
