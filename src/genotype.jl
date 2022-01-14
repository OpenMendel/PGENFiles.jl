const flipmap = [0x02, 0x01, 0x00, 0x03]
const onebitmap = [0x00 => 0x01, 0x00 => 0x02, 0x00 => 0x03, 
    0xff => 0xff, 0x01 => 0x02, 0x01 => 0x03, 0xff => 0xff, 0xff => 0xff, 0x02 => 0x03]

"""
    get_genotypes!(buf, p, v; ldbuf)

Computes unphased biallelic dosage of ALT allele. 

- `buf`: stores genotype values.
- `p`: a `Pgen` object.
- `v`: a `Variant` object.
- `ldbuf`: most recent non-LD-compressed genotypes.

Returns: 
- `buf`: resulting genotype values
- `variant_record`: current record for variant
- `offset`: end of dosage record on the current variant record track.
"""
function get_genotypes!(buf::Vector{UInt8}, p::Pgen, v::Variant; 
    ldbuf::Union{Nothing, Vector{UInt8}}=nothing)
    compression_type = v.record_type & 0x07
    @assert v.record_type & 0x08 == 0 "Multiallelic case unsupported"
    if (compression_type == 0x02 || compression_type == 0x03) 
        if ldbuf === nothing
            # load most recent non-LD-compressed dosage
            get_genotypes!(buf, p, p.header.most_recent_non_ld[v.index])
        else
            buf .= ldbuf
        end
    end
    if p.variant_record_cache !== nothing
        seek(p.io, v.offset)
        p.variant_record_cache .= read(p.io, v.length)
        variant_record = p.variant_record_cache
    else
        variant_record = @view p.data[v.offset + 1 : v.offset + v.length]
    end
    if compression_type == 0x00
        offset = _get_genotypes_no_compression!(buf, p, variant_record)
    elseif compression_type == 0x01
        offset = _get_genotypes_1bit!(buf, p, variant_record)
    elseif compression_type == 0x02 || compression_type == 0x03
        offset = _get_genotypes_difflist!(buf, p, variant_record)
        if compression_type == 0x03
            for i in 1:p.header.n_samples
                buf[i] = flipmap[buf[i] + 0x01]
            end
        end
    elseif compression_type == 0x04 || compression_type == 0x06 || compression_type == 0x07
        if compression_type == 0x04
            fill!(buf, 0x00)
        elseif compression_type == 0x06
            fill!(buf, 0x02)
        else
            fill!(buf, 0x03)
        end
        offset = _get_genotypes_difflist!(buf, p, variant_record)
    else
        @error "invalid compression type"
    end
    buf, variant_record, offset
end

"""
    get_genotypes(p, v)

Computes unphased biallelic dosage of ALT allele. 

- `p`: a `Pgen` object.
- `v`: a `Variant` object.

Returns: 
- `buf`: resulting genotype values
- `variant_record`: current record for variant
- `offset`: end of dosage record on the current variant record track.
"""
function get_genotypes(p::Pgen, v::Variant)
    buf = Vector{UInt8}(undef, p.header.n_samples)
    get_genotypes!(buf, p, v)
end

# for compression mode "0x00", no compression
function _get_genotypes_no_compression!(buf::Vector{UInt8}, p::Pgen, variant_record::AbstractVector{UInt8})
    n_samples = p.header.n_samples
    n_bytes = (n_samples + 3) >> 2
    genotypes_raw_cache = @view variant_record[1:n_bytes]
    @inbounds for i in 1:n_samples
        ip3 = i + 3
        buf[i] = (genotypes_raw_cache[ip3 >> 2] >> ((ip3 & 0x03) << 1)) & 0x03
    end
    n_bytes
end

function _get_genotypes_1bit!(buf::Vector{UInt8}, p::Pgen, variant_record::AbstractVector{UInt8})
    n_samples = p.header.n_samples
    n_bytes = (n_samples + 7) >> 3 
    falseval, trueval = onebitmap[variant_record[1]]
    bv = BitsVector(@view(variant_record[2:2 + n_bytes - 1]), 1, n_samples)
    @inbounds for i in 1:n_samples
        buf[i] = bv[i] == 0x01 ? trueval : falseval
    end
    dl, offset = parse_difflist(variant_record, UInt(n_bytes + 1), p.header.bytes_per_sample_id, true)
    _get_difflist_genotypes!(buf, p, dl)
    offset
end

function _get_genotypes_difflist!(buf::Vector{UInt8}, p::Pgen, variant_record::AbstractVector{UInt8})
    dl, offset = parse_difflist(variant_record, zero(UInt), p.header.bytes_per_sample_id, true)
    _get_difflist_genotypes!(buf, p, dl)
    offset
end

function _get_difflist_genotypes!(buf::Vector{UInt8}, p::Pgen, dl::DiffList)
    ngroups = (dl.len + 63) ÷ 64
    for gid in 1:ngroups
        parse_difflist_sampleids!(p.difflist_cache, p.difflist_cache_incr, dl, gid)
        for (idx, sampleid) in enumerate(p.difflist_cache)
            totalidx = 64 * (gid - 1) + idx
            if totalidx > dl.len
                break
            end
            buf[sampleid] = dl.genotypes[totalidx]
        end
    end
end
