const basemap_Float32 = [0.0f0, 1.0f0, 2.0f0, NaN32]
const basemap_Float64 = [0.0, 1.0, 2.0, NaN]
const dosage_unit_Float64 = 1/(1 << 14)
const dosage_unit_Float32 = Float32(dosage_unit_Float64)

"""
    alt_allele_dosage!(buf, genobuf, p, v; genoldbuf)

Computes unphased biallelic dosage of ALT allele. 

- `buf`: stores dosage values.
- `genobuf`: stores genotype values.
- `p`: a `Pgen` object.
- `v`: a `Variant` object.
- `genoldbuf`: most recent non-LD-compressed genotypes.

Returns: 
- `buf`
- `genobuf`
- `offset`: end of dosage record on the current variant record track.
"""
function alt_allele_dosage!(buf::AbstractVector{T}, genobuf::AbstractVector{UInt8}, 
    p::Pgen, v::Variant;
    genoldbuf::Union{Nothing, Vector{UInt8}}=nothing
    ) where T <: AbstractFloat
    n_samples = p.header.n_samples
    @assert length(buf) == n_samples && length(genobuf) == n_samples

    _, variant_record, offset = get_genotypes!(genobuf, p, v; prev_buf=genoldbuf)

    # convert genotype to floating-point values
    basemap = T == Float32 ? basemap_Float32 : basemap_Float64
    @inbounds for i in 1:n_samples
        buf[i] = basemap[genobuf[i] + 1]
    end

    @assert v.record_type & 0x08 == 0 "Multiallelic case unsupported"
    if v.record_type & 0x08 != 0
        # TODO: compute offset for multiallelic genotype track
    end

    # TODO: test this.
    if v.record_type & 0x10 != 0 # skip phase info (aux track 2) 
        offset = _phase_offset(p, genobuf, variant_record, offset)
    end

    if v.record_type & 0x60 == 0
        return buf, variant_record, offset
    elseif v.record_type & 0x20 !== 0 && v.record_type & 0x40 == 0 
        # track 3 is a difflist
        dl, offset = parse_difflist(variant_record, offset, p.header.bytes_per_sample_id, false)
        offset = _get_difflist_dosage!(buf, p, dl, variant_record, offset)
    elseif v.record_type & 0x20 == 0 && v.record_type & 0x40 != 0 
        # track 3 does not exist, 0xffff is missing value in track 4
        dosages = reinterpret(UInt16, @view(variant_record[offset + 1 : offset + 2 * n_samples]))
        dosage_unit = T == Float32 ? dosage_unit_Float32 : dosage_unit_Float64
        for i in 1:n_samples
            d = dosages[i]
            buf[i] = d == 0xffff ? NaN : d * dosage_unit
        end
        offset += 2 * n_samples
    else # track 3 is a bitarray for dosage existence in track 4. 
        bv_nbytes = (n_samples + 7) >> 3
        bv = BitsVector(@view(variant_record[offset + 1 : offset + bv_nbytes]), 1, n_samples)
        offset += bv_nbytes
        dosage_len = sum(bv)
        dosages = reinterpret(UInt16, @view(variant_record[offset + 1 : offset + 2 * dosage_len]))
        dosage_unit = T == Float32 ? dosage_unit_Float32 : dosage_unit_Float64
        j = 1
        for i in 1:n_samples
            if bv[i] == 0x01
                buf[i] = dosages[j] * dosage_unit
                j += 1
            end
            if j > dosage_len
                break
            end
        end
        offset += 2 * dosage_len
    end
    return buf, genobuf, offset
end

"""
    ref_allele_dosage!(buf, genobuf, p, v; genoldbuf)

Computes unphased biallelic dosage of REF allele. 

- `buf`: stores dosage values.
- `genobuf`: stores genotype values.
- `p`: a `Pgen` object.
- `v`: a `Variant` object.
- `genoldbuf`: most recent non-LD-compressed genotypes.

Returns: 
- `buf`
- `genobuf`
- `offset`: end of dosage record on the current variant record track.
"""
function ref_allele_dosage!(buf::AbstractVector{T}, genobuf::AbstractVector{UInt8}, 
    p::Pgen, v::Variant; 
    genoldbuf::Union{Nothing, Vector{UInt8}}=nothing
    ) where T <: AbstractFloat
    _, offset = alt_allele_dosage!(buf, genobuf, p, v; genoldbuf=genoldbuf)
    @inbounds for i in 1:p.header.n_samples
        buf[i] = 2.0 - buf[i]
    end
    return buf, genobuf, offset
end

"""
    alt_allele_dosage(p, v)

Computes unphased biallelic dosage of ALT allele. 

- `p`: a `Pgen` object.
- `v`: a `Variant` object.

Returns: 
- `buf`: stores dosage values.
- `genobuf`: stores genotype values.
- `offset`: end of dosage record on the current variant record track.
"""
function alt_allele_dosage(p::Pgen, v::Variant)
    n_samples = p.header.n_samples
    buf = Vector{Float32}(undef, n_samples)
    genobuf = Vector{UInt8}(undef, n_samples)
    alt_allele_dosage!(buf, genobuf, p, v)
end

"""
    ref_allele_dosage(p, v)

Computes unphased biallelic dosage of REF allele. 

- `p`: a `Pgen` object.
- `v`: a `Variant` object.

Returns: 
- `buf`: stores dosage values.
- `genobuf`: stores genotype values.
- `offset`: end of dosage record on the current variant record track.
"""
function ref_allele_dosage(p::Pgen, v::Variant)
    n_samples = p.header.n_samples
    buf = Vector{Float32}(undef, n_samples)
    genobuf = Vector{UInt8}(undef, n_samples)
    ref_allele_dosage!(buf, genobuf, p, v)
end

# This function is untested.
function _phase_offset(p::Pgen, g::AbstractVector{UInt8},
    variant_record::AbstractVector{UInt8}, offset::Integer)
    hetero_count = 0
    @inbounds for i in 1:p.header.n_samples
        if g[i] == 0x01
            hetero_count += 1
        end
    end
    phasepresent = variant_record[offset + 1] & 0x01 != 0
    phase_bits = 0
    if phasepresent
        phasepresentbv_nbytes = (hetero_count >> 3) + 1 #((hetero_count + 1) + 7) >> 3
        phase_bits = 0
        for i in 1:phasepresentbv_nbytes
            phase_bits += count_ones(variant_record[offset + i])
        end
        phase_bits -= 1 # remove initial set bit.
        offset += phasepresentbv_nbytes 
    else
        phase_bits = hetero_count + 1 # should include phasepresent bit
    end
    phase_nbytes = (phase_bits + 7) >> 3
    offset += phase_nbytes
    return offset
end

# replaces the dosage in `buf` with the values in the `DiffList` `dl`. 
function _get_difflist_dosage!(buf::Vector{T}, p::Pgen, dl::DiffList, 
    variant_record::AbstractVector{UInt8}, offset::Integer) where T <: AbstractFloat
    ngroups = (dl.len + 63) ÷ 64
    if dl.len == 0 
        return offset
    end
    dosage_unit = T == Float32 ? dosage_unit_Float32 : dosage_unit_Float64
    dosages = reinterpret(UInt16, @view(variant_record[offset + 1 : offset + 2 * dl.len]))
    for gid in 1:ngroups
        parse_difflist_sampleids!(p.difflist_cache, p.difflist_cache_incr, dl, gid)
        for (idx, sampleid) in enumerate(p.difflist_cache)
            totalidx = 64 * (gid - 1) + idx
            if totalidx > dl.len
                break
            end
            buf[sampleid] = dosages[totalidx] * dosage_unit
        end
    end
    return offset + 2 * dl.len
end
