struct Header
    # magic number (first two bytes): 0x6c 0x1b. 
    storage_mode::UInt # third byte. restrict to 0x10 for now.
    n_variants::UInt # 4th-7th byte. 
    n_samples::UInt # 8th-11th byte. 
    bits_per_variant_type::UInt # bits 0-3 of 12th byte
    bytes_per_record_length::UInt # bits 0-3 of 12th byte
    bytes_per_allele_count::UInt # bits 4-5 of 12th byte, restrict to 0 for now (no multiallelic variants allowed).
    provisional_reference::UInt # bits 6-7 of 12th byte. 
    n_blocks::UInt # number of blocks of 2^16 variants. Int(ceil(n_variants / 2 ^ 16)). 
    variant_block_offsets::Vector{UInt} # record starting points of #0, #65536, ... length of (8 * n_blocks) bytes.
    # The following appear in blocks of 2^16 variants.
    variant_types::Union{ScatteredBitsVector, ScatteredVector}
    variant_lengths::ScatteredVector
    allele_counts::Union{ScatteredVector, Nothing}
    provisional_reference_flags::Union{ScatteredBitsVector, Nothing}
end

struct Pgen
    io::IOStream
    data::Vector{UInt8}
    header::Header
    genotypes_cache::BitsVector{Vector{UInt8}} # 2-bit packed array of length header.n_samples
    difflist_cache::Vector{UInt8} # length-64 vector
end
