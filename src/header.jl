function Header(data::Vector{UInt8})
    #seek(io, 0)
    # check magic number
    (data[1] == 0x6c && data[2] == 0x1b) || throw(ArgumentError("wrong magic number in the input file"))
    storage_mode = data[3]
    storage_mode == 0x10 || throw(ArgumentError("Only the standard format supported for now"))
    n_variants = reinterpret(UInt32, view(data, 4:7))[1]
    n_samples = reinterpret(UInt32, view(data, 8:11))[1]
    format_num = data[12]

    # figure out bits_per_variant_type and bytes_per_record_length
    format_num_0_3 = format_num & 0x0f

    (bits_per_variant_type, bytes_per_record_length) = variant_type_lengths[format_num_0_3]

    bytes_allele_counts = (format_num & 0x30) >> 4

    # provisional reference
    provisional_reference = (format_num & 0xc0) >> 6

    n_blocks = ceil_int(n_variants, 2 ^ 16) #Int(ceil(n_variants / 2 ^ 16))

    variant_block_offsets = reinterpret(UInt64, view(data, 13:(12 + n_blocks * 8)))

    #remaining_header = mmap(io, Vector{UInt8}, variant_block_offsets[1] - position(io))
    offset = convert(UInt64, 12 + n_blocks * 8)

    t_variant_types = nothing
    sectors_variant_types = []

    dtype_variant_length = bytes_to_UInt[bytes_per_record_length]
    t_variant_sizes = nothing
    sectors_variant_lengths = []

    if bytes_allele_counts != 0
        dtype_allele_counts = bytes_to_UInt[bytes_allele_counts]
        t_allele_counts = nothing
        sectors_allele_counts = []
    else
        sectors_allele_counts = nothing
    end
    if provisional_reference == 0x03
        t_provisional_reference_flags = nothing
        sectors_provisional_reference = []
    else
        sectors_provisional_reference = nothing
    end
    for i in 1:n_blocks
        block_size = i < n_blocks ? convert(UInt, 2 ^ 16) : begin
            rem = n_variants % convert(UInt, 2 ^ 16)
            rem > 0 ? rem : convert(UInt, 2 ^ 16)
        end

        # read variant type track
        data_per_byte = 8 ÷ bits_per_variant_type
        size_variant_types = ceil_int(block_size, data_per_byte)
        #convert(UInt, ceil(block_size / (8 ÷ bits_per_variant_type)))
        arr = view(data, 
            offset + 1 : offset + size_variant_types)
        push!(sectors_variant_types, Ref(arr))
        t_variant_types = typeof(arr)
        offset += size_variant_types

        # read variant record length track
        arr = view(
            data, offset + 1 : 
            offset + block_size * bytes_per_record_length)
        reinterpreted = reinterpret(dtype_variant_length, arr)
        t_variant_sizes = typeof(reinterpreted)
        push!(sectors_variant_lengths, Ref(reinterpreted))
        offset += bytes_per_record_length * block_size 

        # read allele counts track
        if sectors_allele_counts !== nothing
            arr = view(
                data, offset + 1:
                offset + block_size * bytes_allele_counts)
            reinterpreted = reinterpret(dtype_allele_counts, arr)
            t_allele_counts = typeof(reinterpreted)
            push!(sectors_allele_counts, Ref(reinterpreted))
            offset += bytes_allele_counts * block_size
        end

        # read provisional reference flags track
        if sectors_provisional_reference !== nothing
            size_pr = ceil_int(block_size, 8)#convert(UInt, ceil(block_size / 8))
            arr = view(
                data, offset + 1 :
                offset + size_pr)
            t_provisional_reference_flags = typeof(arr)
            push!(sectors_provisional_reference, Ref(arr))
            offset += size_pr
        end
    end
    if bits_per_variant_type == 8
        variant_types = ScatteredVector{UInt8, t_variant_types}(
            sectors_variant_types, n_variants, n_blocks)
    else
        variant_types = ScatteredBitsVector{t_variant_types}(
            sectors_variant_types, bits_per_variant_type, 
            n_variants, n_blocks)
    end
    variant_lengths = ScatteredVector{dtype_variant_length,t_variant_sizes}(
        sectors_variant_lengths, 
        n_variants, n_blocks)
    if sectors_allele_counts !== nothing
        allele_counts = ScatteredVector{dtype_allele_counts,t_allele_counts}(
            sectors_allele_counts, 
            n_variants, n_blocks)
    else
        allele_counts = nothing
    end
    if sectors_provisional_reference !== nothing
        provisional_reference_flags = ScatteredBitsVector{t_provisional_reference_flags}(
            sectors_provisional_reference, 1, 
            n_variants, n_blocks)
    else
        provisional_reference_flags = nothing
    end
    Header(storage_mode, n_variants, n_samples, 
        bits_per_variant_type, bytes_per_record_length, 
        bytes_allele_counts, provisional_reference, 
        n_blocks, variant_block_offsets, 
        variant_types, variant_lengths, allele_counts, provisional_reference_flags)
end

function Header(filename::String)
    io = open(filename)
    data = mmap(io)
    Header(data)
end
