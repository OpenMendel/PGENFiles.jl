#TODO: test it.
function DiffList(data::AbstractVector{UInt8}, 
    offset::Base.RefValue{UInt}, 
    bytes_per_sampleid::Integer,
    has_genotype::Bool)
    # length
    len = decode_single(data; offset=offset)
    println(len)
    # sample id bases
    n_groups = ceil_int(len, UInt32(64))
    sample_id_dtype = bytes_to_UInt[bytes_per_sampleid]
    sample_id_bases = reinterpret(sample_id_dtype, view(data, 
        (offset[] + 1): (offset[] + n_groups * bytes_per_sampleid)))
    offset[] += n_groups * bytes_per_sampleid

    # sizes of the final components
    final_component_sizes = view(data, offset[] + 1: offset[] + n_groups - 1)
    offset[] += n_groups - 1

    # genotypes
    if has_genotype
        genotype_bytes = ceil_int(len, 4)
        genotype_data = view(data, offset[] + 1 : offset[] + genotype_bytes)
        genotypes = BitsVector{typeof(genotype_data)}(Ref(genotype_data), 2, len)
        offset[] += genotype_bytes
    else
        genotypes = nothing
    end

    # final component: differences of indices
    final_component_size = sum(final_component_sizes) + 63 * length(final_component_sizes)
    last_group_size = len % 64 - 1
    final_component_size += size_n(data, last_group_size, offset[] + final_component_size)
    final_component = view(data, offset[] + 1 : offset[] + final_component_size)
    offset[] += final_component_size

    DiffList{typeof(sample_id_bases), typeof(final_component_sizes), 
        typeof(genotype_data), typeof(final_component)}(
        len, Ref(sample_id_bases), Ref(final_component_sizes), 
        has_genotype,
        genotypes, Ref(final_component)), offset[]
end