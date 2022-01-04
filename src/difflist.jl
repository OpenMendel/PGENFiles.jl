function parse_difflist(data::AbstractVector{UInt8}, 
    offset::UInt, 
    bytes_per_sample_id::Integer,
    has_genotype::Bool)
    # length
    len, offset = decode_single(data; offset=offset)

    # sample id bases
    n_groups = ceil_int(len, 0x000040) # 64 in decimal
    sample_id_dtype = bytes_to_UInt[bytes_per_sample_id]
    sample_id_bases = reinterpret(sample_id_dtype, view(data, 
        (offset + 1): (offset + n_groups * bytes_per_sample_id)))
    offset += n_groups * bytes_per_sample_id

    # sizes of the final components
    final_component_sizes = view(data, offset + 1: offset + n_groups - 1)
    offset += n_groups - 1

    # genotypes
    if has_genotype
        genotype_bytes = ceil_int(len, 4)
        genotype_data = view(data, offset + 1 : offset + genotype_bytes)
        genotypes = BitsVector{typeof(genotype_data)}(Ref(genotype_data), 2, len)
        offset += genotype_bytes
    else
        genotypes = nothing
    end

    # final component: differences of indices
    final_component_size = sum(final_component_sizes) + 63 * length(final_component_sizes)
    last_group_size = len % 64 - 1
    last_group_size = last_group_size == -1 ? 63 : last_group_size
    last_incr_size = size_n(data, last_group_size, offset + final_component_size)
    final_component_size += last_incr_size
    final_component = view(data, offset + 1 : offset + final_component_size)
    offset += final_component_size

    DiffList{typeof(sample_id_bases), typeof(final_component_sizes),
        typeof(genotypes), typeof(final_component)}(
        len, Ref(sample_id_bases), Ref(final_component_sizes), 
        has_genotype,
        genotypes, Ref(final_component)), offset
end

function parse_difflist_sampleids!(idx::AbstractArray, idx_incr::AbstractArray, 
        dl::DiffList, gid::Integer, sid_incr_offset::Union{Nothing, UInt} = nothing)
    n_groups = ceil_int(dl.len, 64)
    @assert gid <= n_groups
    if sid_incr_offset === nothing
        sid_incr_offset = zero(UInt)
        @inbounds for i in 1:(gid - 1)
            sid_incr_offset += dl.last_component_sizes[][i] + 0x3f # offset by 63 (0x3f)
        end
    end
    idx_incr[1] = 0
    baseidx = dl.sample_id_bases[][gid]
    if gid != n_groups
        decode_multiple!(@view(idx_incr[2:end]), dl.sample_id_increments[]; offset = UInt(sid_incr_offset))
        cumsum!(idx_incr, idx_incr)
        idx .= baseidx .+ idx_incr
    else
        count = dl.len % 64
        count == 0 ? 64 : count
        decode_multiple!(@view(idx_incr[2:end]), dl.sample_id_increments[]; 
            count = count - 1, offset = UInt(sid_incr_offset))
        cumsum!(view(idx_incr, 1:count), view(idx_incr, 1:count))
        idx[1:count] .= baseidx .+ view(idx_incr, 1:count)
        fill!(@view(idx[(count + 1):end]), 0)
    end
    idx
end
