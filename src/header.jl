function Header(io::IOStream)
    seek(io, 0)
    # check magic number
    read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in the input file"))
    storage_mode = read(io, UInt8)
    storage_mode == 0x10 || throw(ArgumentError("Only the standard format supported for now"))
    n_variants = Int(read(io, UInt32))
    n_samples = Int(read(io, UInt32))
    format_num = read(io, UInt8)

    # figure out bits_per_variant_type and bytes_per_record_length
    format_num_0_3 = format_num & 0x0f
    variant_type_lengths = Dict(
        0x00 => (4, 1), 0x01 => (4, 2), 0x02 => (4, 3), 0x03 => (4, 4),
        0x04 => (8, 1), 0x05 => (8, 2), 0x06 => (8, 3), 0x07 => (8, 4)
    )
    (bits_per_variant_type, bytes_per_record_length) = variant_type_lengths[format_num_0_3]

    bytes_allele_counts = (format_num & 0x30) >> 4

    # provisional reference
    provisional_reference = (format_num & 0xc0) >> 6

    n_blocks = Int(ceil(n_variants / 2 ^ 16))

    variant_block_offsets = read(io, n_blocks * sizeof(UInt64))
    variant_block_offsets = reinterpret(UInt64, variant_starts)

    sectors_variant_types = Vector{UInt8}[]
    dtype_variant_length = bytes_to_UInt[bytes_per_record_length]
    sectors_variant_lengths = Vector{dtype_variant_length}[]

    if bytes_allele_counts != 0
        dtype_allele_counts = bytes_to_UInt[bytes_allele_counts]
        sectors_allele_counts = Vector{dtype_allele_counts}[]
    else
        sectors_allele_counts = nothing
    end
    if provisional_reference != 0x03
        sectors_provisional_reference = Vector{UInt8}[]
    else
        sectors_provisional_reference = nothing
    end
    offset = position(io)
    for i in 1:n_blocks
        block_size = i < n_blocks ? 2 ^ 16 : begin
            rem = n_variants % 2 ^ 16
            rem > 0 ? rem : 2 ^ 16
        end
        size_variant_types = convert(Int, ceil(block_size / (8 ÷ bits_per_variant_type)))
        push!(mmap(sectors_variant_types, Vector{UInt8}, (size_variant_types,), offset))
        offset += size_variant_types
        push!(mmap(sectors_variant_lengths, Vector{dtype_variant_length}, 
            (block_size,), offset))
        offset += bytes_per_record_length * block_size 
        if sectors_allele_counts !== nothing
            push!(mmap(sectors_allele_counts, Vector{dtype_allele_counts}, 
                (block_size,), offset))
            offset += bytes_allele_counts * block_size
        end
        if sectors_provisional_reference !== nothing
            size_pr = convert(Int, ceil(block_size / 8))
            push!(mmap(sectors_provisional_reference, Vector{UInt8}, 
                (size_pr,), offset))
            offset += size_pr
        end
    end
    if bits_per_variant_type == 8
        variant_types = ScatteredVector{UInt8}(sectors_variant_types, n_variants, n_blocks)
    else
        variant_types = ScatteredBitsVector(sectors_variant_types, bits_per_variant_type, 
            n_variants, n_blocks)
    end
    variant_lengths = ScatteredVector{dtype_variant_length}(sectors_variant_lengths, 
        n_variants, n_blocks)
    if sectors_allele_counts !== nothing
        allele_counts = ScatteredVector{dtype_allele_counts}(sectors_allele_counts, 
            n_variants, n_blocks)
    else
        allele_counts = nothing
    end
    if sectors_provisional_reference !== nothing
        provisional_reference_flags = ScatteredBitsVector(sectors_provisional_reference, 1, 
            n_variants, n_blocks)
    else
        provisional_reference_flags = nothing
    end

    Header(storage_mode, n_variants, n_samples, 
        bits_per_variant_type, bytes_per_record_length, 
        bytes_per_allele_count, provisional_reference, 
        n_blocks, variant_block_offsets, 
        variant_types, variant_lengths, allele_counts, provisional_reference_flags)
end

function Header(filename::String)
    io = open(filename)
    Header(io)
end
