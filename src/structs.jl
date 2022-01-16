mutable struct Variant
    index::UInt64 # 1-based
    offset::UInt64
    record_type::UInt8
    length::UInt64
end

struct Header{VTT,VLT,ACT,PRFT}
    # magic number (first two bytes): 0x6c 0x1b. 
    storage_mode::UInt # third byte. restrict to 0x10 for now.
    n_variants::UInt # 4th-7th byte. 
    n_samples::UInt # 8th-11th byte. 
    bits_per_variant_type::UInt # bits 0-3 of 12th byte
    bytes_per_record_length::UInt # bits 0-3 of 12th byte
    bytes_per_allele_count::UInt # bits 4-5 of 12th byte, restrict to 0 for now (no multiallelic variants allowed).
    bytes_per_sample_id::UInt # determined by n_samples (1 if < 2 ^ 8, etc.)
    provisional_reference::UInt # bits 6-7 of 12th byte. 
    n_blocks::UInt # number of blocks of 2^16 variants. Int(ceil(n_variants / 2 ^ 16)). 
    variant_block_offsets::Vector{UInt} # record starting points of #0, #65536, ... length of (8 * n_blocks) bytes.
    # The following appear in blocks of 2^16 variants.
    variant_types::VTT # Union{ScatteredBitsVector, ScatteredVector}
    variant_lengths::VLT # ScatteredVector
    allele_counts::ACT # Union{ScatteredVector, Nothing}
    provisional_reference_flags::PRFT # Union{ScatteredBitsVector, Nothing}
    most_recent_non_ld::Dict{UInt, Variant}
end

struct Pgen{ST}
    io::IOStream
    data::Union{Nothing, Vector{UInt8}}
    header::Header
    variant_record_cache::Union{Nothing, Vector{UInt8}} # used only with no_mmap
    difflist_cache::Vector{ST} # length-64 vector for 64 Sample IDs.
    difflist_cache_incr::Vector{UInt32}
end

"""
    Pgen(filename; no_mmap)

Creates an instance of `Pgen` from `filename`. `no_mmap` chooses whether to use `Mmap`.
"""
function Pgen(filename::String; no_mmap::Bool=false)
    io = open(filename)
    if !no_mmap
        data = mmap(io)
    else
        data = nothing
    end
    header = Header(io)
    ST = bytes_to_UInt[header.bytes_per_sample_id]
    if !no_mmap
        variant_record_cache = nothing
    else
        variant_record_cache = Vector{UInt8}(undef, maximum(header.variant_lengths))
    end
    difflist_cache = Vector{ST}(undef, 64)
    difflist_cache_incr = Vector{UInt32}(undef, 64)
    difflist_cache_incr[1] = 0
    Pgen{ST}(io, data, header, variant_record_cache,
        difflist_cache, difflist_cache_incr)
end

@inline n_variants(p::Pgen) = p.header.n_variants
@inline n_samples(p::Pgen) = p.header.n_samples 
