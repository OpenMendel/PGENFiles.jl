"""
    parse_difflist(data, offset, bytes_per_sample_id, has_genotypes)

Parses a single `DiffList` from `data`, starting with `offset + 1`-th byte.
"""
function parse_difflist(data::AbstractVector{UInt8}, 
    offset::Integer, 
    bytes_per_sample_id::Integer,
    has_genotype::Bool)
    # length
    len, offset = decode_single(pointer(data); offset=offset)
    # length-zero list
    sample_id_dtype = bytes_to_UInt[bytes_per_sample_id]
    if len == 0
        return DiffList{sample_id_dtype, Nothing}(
            0, C_NULL, C_NULL, has_genotype, nothing, C_NULL), offset
    end
    # sample id bases
    n_groups = ceil_int(len, 0x000040) # 64 in decimal
    sample_id_bases = reinterpret(sample_id_dtype, view(data, 
        Int(offset + 1) : Int(offset + n_groups * bytes_per_sample_id)))

    offset += n_groups * bytes_per_sample_id

    # sizes of the final components
    final_component_sizes = view(data, Int(offset + 1): Int(offset + n_groups - 1))
    offset += n_groups - 1

    # genotypes
    if has_genotype
        genotype_bytes = ceil_int(len, 4)
        # genotype_data = view(data, offset + 1 : offset + genotype_bytes)
        genotypes = BitsVector(pointer(data, offset + 1), 2, len)
        # BitsVector{typeof(genotype_data)}(Ref(genotype_data), 2, len)
        offset += genotype_bytes
    else
        genotypes = nothing
    end

    # final component: differences of indices
    final_component_size = sum(final_component_sizes) + 63 * length(final_component_sizes)
    last_group_size = len % 64
    last_group_size = last_group_size == 0 ? 64 : last_group_size
    # increment has one less value
    last_group_size -= 1
    last_incr_size = size_n(data, last_group_size, offset + final_component_size)
    final_component_size += last_incr_size
    final_component = pointer(data, offset + 1)#view(data, Int(offset + 1) : Int(offset + final_component_size))

    # TODO: defer computation of offset of the last block to parse_difflist_sampleids!
    offset += final_component_size

    DiffList{typeof(sample_id_bases),typeof(genotypes)}(
        len, sample_id_bases, pointer(final_component_sizes), 
        has_genotype,
        genotypes, final_component), offset
end

"""
    parse_difflist_sampleids!(idx, idx_incr, dl, gid, sid_incr_offset=nothing)

Parses a `gid`-th group of 64 sample ids in a `DiffList` `dl` into `idx`. 
"""
function parse_difflist_sampleids!(idx::AbstractArray, idx_incr::AbstractArray, 
        dl::DiffList, gid::Integer, sid_incr_offset::Union{Nothing, UInt} = nothing)
    # TODO: compute offset of the last block in this function.
    n_groups = ceil_int(dl.len, 64)
    @assert gid <= n_groups
    if sid_incr_offset === nothing
        sid_incr_offset = zero(UInt)
        @inbounds for i in 1:(gid - 1)
            sid_incr_offset += unsafe_load(dl.last_component_sizes, i) + 0x3f
            #dl.last_component_sizes[][i] + 0x3f # offset by 63 (0x3f)
        end
    end
    idx_incr[1] = 0
    # offset by one to make it 1-based.
    baseidx = dl.sample_id_bases[gid] + 1
    if gid != n_groups
        decode_multiple!(pointer(idx_incr, 2), dl.sample_id_increments; 
            count = 63, offset = UInt(sid_incr_offset))
        cumsum!(idx_incr, idx_incr)
        idx .= baseidx .+ idx_incr
    else
        count = dl.len % 64
        count = count == 0 ? 0x00000040 : count
        decode_multiple!(pointer(idx_incr, 2), dl.sample_id_increments; 
            count = count - 1, offset = UInt(sid_incr_offset))
        idx[1] = baseidx
        @inbounds for i in 2:count
            idx[i] = idx[i-1] + idx_incr[i]
        end
        @inbounds for i in (count + 1):length(idx)
            idx[i] = 0
        end
        # cumsum!(view(idx_incr, 1:count), view(idx_incr, 1:count))
        # idx[1:count] .= baseidx .+ view(idx_incr, 1:count)
        # fill!(@view(idx[(count + 1):end]), 0)
    end
    idx
end
