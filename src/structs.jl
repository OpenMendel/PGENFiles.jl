struct Header
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
    variant_types::Union{ScatteredBitsVector, ScatteredVector}
    variant_lengths::ScatteredVector
    allele_counts::Union{ScatteredVector, Nothing}
    provisional_reference_flags::Union{ScatteredBitsVector, Nothing}
end

struct Pgen{ST}
    io::IOStream
    data::Vector{UInt8}
    header::Header
    genotypes_prev::Vector{UInt8} # for LD-compressed genotypes
    genotypes_cache::Vector{UInt8} # 0x00, 0x01, or 0x02. Byte-aligned for performance.
    dosage_cache::Vector{Float32} # Dosage values are represented by 16-bit numbers, Float32 is enough. 
    difflist_cache::Vector{ST} # length-64 vector for 64 Sample IDs.
    difflist_cache_incr::Vector{UInt32}
end

function Pgen(filename::String)
    io = open(filename)
    data = mmap(io)
    header = Header(data)
    ST = bytes_to_UInt[header.bytes_per_sample_id]
    genotypes_prev = Vector{UInt8}(undef, header.n_samples)
    genotypes_cache = Vector{UInt8}(undef, header.n_samples)
    dosage_cache = Vector{Float32}(undef, header.n_samples)
    difflist_cache = Vector{ST}(undef, 64)
    difflist_cache_incr = Vector{UInt32}(undef, 64)
    difflist_cache_incr[1] = 0
    Pgen{ST}(io, data, header, genotypes_prev, genotypes_cache, dosage_cache, 
        difflist_cache, difflist_cache_incr)
end
