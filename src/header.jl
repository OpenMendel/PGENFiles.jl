function Header(io::IOStream)
    # check magic number
    (read(io, UInt8) == 0x6c && read(io, UInt8) == 0x1b) || throw(ArgumentError("wrong magic number in the input file"))
    # 2 bytes so far
    storage_mode = read(io, UInt8) # 3 bytes
    storage_mode == 0x10 || throw(ArgumentError("Only the standard format supported for now"))
    n_variants = read(io, UInt32) # 7 bytes
    n_samples = read(io, UInt32) # 11 bytes

    if n_samples <= 2 ^ 8
        bytes_per_sample_id = 1
    elseif n_samples <= 2 ^ 16
        bytes_per_sample_id = 2
    elseif n_samples <= 2 ^ 24
        bytes_per_sample_id = 3
    else
        bytes_per_sample_id = 4
    end
    
    format_num = read(io, UInt8) # 12 bytes

    # figure out bits_per_variant_type and bytes_per_record_length
    format_num_0_3 = format_num & 0x0f

    (bits_per_variant_type, bytes_per_record_length) = variant_type_lengths[format_num_0_3]

    bytes_allele_counts = (format_num & 0x30) >> 4

    # provisional reference
    provisional_reference = (format_num & 0xc0) >> 6

    n_blocks = ceil_int(n_variants, 2 ^ 16) #Int(ceil(n_variants / 2 ^ 16))
    variant_block_offsets = Vector{UInt64}(undef, n_blocks)
    read!(io, variant_block_offsets) # 12 + 8n_blocks
    #variant_block_offsets = reinterpret(UInt64, data[13:(12 + n_blocks * 8)])

    offset = convert(UInt64, 12 + 8n_blocks)

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
        arr = read(io, size_variant_types) #data[offset + 1 : offset + size_variant_types]
        push!(sectors_variant_types, Ref(arr))
        t_variant_types = typeof(arr)
        offset += size_variant_types

        # read variant record length track
        arr = read(io, block_size * bytes_per_record_length) #data[offset + 1 : offset + block_size * bytes_per_record_length]
        reinterpreted = reinterpret(dtype_variant_length, arr)
        t_variant_sizes = typeof(reinterpreted)
        push!(sectors_variant_lengths, Ref(reinterpreted))
        offset += bytes_per_record_length * block_size 

        # read allele counts track
        if sectors_allele_counts !== nothing
            arr = read(io, block_size * bytes_allele_counts)
            #arr = data[offset + 1:offset + block_size * bytes_allele_counts]
            reinterpreted = reinterpret(dtype_allele_counts, arr)
            t_allele_counts = typeof(reinterpreted)
            push!(sectors_allele_counts, Ref(reinterpreted))
            offset += bytes_allele_counts * block_size
        end

        # read provisional reference flags track
        if sectors_provisional_reference !== nothing
            size_pr = ceil_int(block_size, 8)#convert(UInt, ceil(block_size / 8))
            arr = read(io, size_pr)
            #arr = data[offset + 1 : offset + size_pr]
            t_provisional_reference_flags = typeof(arr)
            push!(sectors_provisional_reference, Ref(arr))
            offset += size_pr
        end
    end
    if bits_per_variant_type == 8
        VTT = ScatteredVector{UInt8, t_variant_types}
        variant_types = ScatteredVector{UInt8, t_variant_types}(
            sectors_variant_types, n_variants, n_blocks)
    else
        VTT = ScatteredBitsVector{t_variant_types}
        variant_types = ScatteredBitsVector{t_variant_types}(
            sectors_variant_types, bits_per_variant_type, 
            n_variants, n_blocks)
    end
    VLT = ScatteredVector{dtype_variant_length,t_variant_sizes}
    variant_lengths = ScatteredVector{dtype_variant_length,t_variant_sizes}(
        sectors_variant_lengths, 
        n_variants, n_blocks)
    if sectors_allele_counts !== nothing
        ACT = ScatteredVector{dtype_allele_counts,t_allele_counts}
        allele_counts = ScatteredVector{dtype_allele_counts,t_allele_counts}(
            sectors_allele_counts, 
            n_variants, n_blocks)
    else
        ACT = Nothing
        allele_counts = nothing
    end
    if sectors_provisional_reference !== nothing
        PRFT = ScatteredBitsVector{t_provisional_reference_flags}
        provisional_reference_flags = ScatteredBitsVector{t_provisional_reference_flags}(
            sectors_provisional_reference, 1, 
            n_variants, n_blocks)
    else
        PRFT = Nothing
        provisional_reference_flags = nothing
    end

    # Store the most recent non-LD-compressed variant for each LD-compressed variant
    # Useful for random-access of variants
    most_recent_non_ld = Dict{UInt, Variant}()
    buf = nothing
    offset = variant_block_offsets[1]
    for j in 1:n_variants
        vt = variant_types[j]
        vl = variant_lengths[j]
        if !((vt & 0x07 == 0x02) || (vt & 0x07 == 0x03))
            buf = Variant(j, offset, vt, vl)
        else
            most_recent_non_ld[j] = buf
        end   
        offset += vl
    end

    Header{VTT, VLT, ACT, PRFT}(storage_mode, n_variants, n_samples, 
        bits_per_variant_type, bytes_per_record_length, 
        bytes_allele_counts, bytes_per_sample_id, provisional_reference, 
        n_blocks, variant_block_offsets, 
        variant_types, variant_lengths, allele_counts, provisional_reference_flags,
        most_recent_non_ld)
end

function Header(filename::String)
    io = open(filename)
    Header(io)
end
