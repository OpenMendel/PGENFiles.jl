"""
    write_PGEN(pgen_filename, x, sampleID, variantID)

Saves numeric matrix `x` into a PGEN formatted file. We assume `x`
stores dosage diploid genotype (i.e. `x[i, j] ∈ [0, 2]`). As such, the code 
below implements a "fixed-width storage mode" which can write a standalone
PGEN file in a single sequential pass. 

# Inputs
+ `pgen_filename`: Output filename of PGEN file
+ `x`: Numeric matrix, each row is a sample and each column is a SNPs. It is assumed
    `x[i, j] ∈ [0, 2]` and missing values are represented either with
    `x[i,j] === missing` or `x[i,j] === NaN`.
+ `sampleID`: Vector storing each sample's ID
+ `variantID`: Vector storing each variant's ID

# PGEN format structure
----Header-----
1. Magic number, storage mode, data dimension
2. 12th byte: structure of variant record types/length, allele counts, provision ref
3. Variant block offsets (starting positions of each variant record)
4. A bunch of blocks, each block storing 2^16 variant record types followed by 2^16 
    variant record length
----Variant records-----
5. Variant record #0 starts here
"""
function write_PGEN(
    pgen_filename::AbstractString, 
    x::AbstractMatrix,
    # sampleID::AbstractVector,
    # variantID::AbstractVector
    )
    n_samples, n_variants = size(x)
    n_blocks = PGENFiles.ceil_int(n_variants, 2^16) #Int(ceil(n_variants / 2 ^ 16))
    bytes_written = 0
    # main PGEN file
    open(pgen_filename, "w") do io
        #
        # Construct header
        #
        # magic number
        bytes_written += write(io, 0x6c, 0x1b)
        # storage mode
        bytes_written += write(io, 0x10) # 2 bytes so far (note: 0-based indexing according to manual) 
        # data dimension
        bytes_written += write(io, UInt32(n_variants)) # 6 bytes
        bytes_written += write(io, UInt32(n_samples)) # 10 bytes
        # 11th byte indicating how data is stored
        bits_per_record_type = 8 # need 8 to encode dosages
        bytes_per_record_length = 2 # each variant stored as UInt16 (i.e. requiring 2 bytes)
        twelfth_byte_bits =  "10" * "00" * "0101" # all ref alleles provisional; no allele counts; 8 bits per record type, 2 bytes per record length; 
        bytes_written += write(io, bitstring2byte(twelfth_byte_bits))
        # some constants for computing variant block offsets (i.e. start position for each block)
        variant_offset = 12 + 8n_blocks
        variant_type_offset = Int(2^16 / (8 / bits_per_record_type))
        variant_length_offset = 2^16 * bytes_per_record_length
        variant_offset += (n_blocks - 1) * (variant_type_offset + variant_length_offset)
        last_block_variants = n_variants % 2^16
        variant_offset += Int(last_block_variants / (8 / bits_per_record_type)) + 
            last_block_variants * bytes_per_record_length # 99325313 if bits_per_record_type = 4; bytes_per_record_length = 2; n_samples = 1092; n_variants = 39728178
        # store variant offsets for each block, assuming each variant record has fixed width
        bytes_per_variant_record = 2^16 * bytes_per_record_length * n_samples
        for b in 1:n_blocks
            block_offset = variant_offset + (b - 1) * bytes_per_variant_record
            for x in int2bytes(block_offset, len=8)
                bytes_written += write(io, x)
            end
        end
        # store variant record types and variant record lengths for each block.
        # Here record type is:
        #   "000" (no compression) +
        #   "0" (no multi allelic hard calls) +
        #   "0" (no phased hetero hard calls) +
        #   "01" (dosage exists for all samples, value of 65535 represents missing) +
        #   "0" (explicit phased-dosages abscent)
        variant_record_type = bitstring2byte("00000010")
        variant_record_byte_length = int2bytes(bytes_per_record_length * n_samples, len=2)
        for b in 1:(n_blocks - 1)
            for snp in 1:2^16
                bytes_written += write(io, variant_record_type)
            end
            for snp in 1:2^16
                bytes_written += write(io, variant_record_byte_length)
            end
        end
        # last block 
        remainders = n_variants % 2^16
        for snp in 1:remainders
            bytes_written += write(io, variant_record_type)
        end
        for snp in 1:remainders
            bytes_written += write(io, variant_record_byte_length)
        end
        #
        # Construct variant records 
        #
        for j in 1:n_variants
            bytes_written += write_variant_record(io, @view(x[:, j]))
        end
    end
    # handle psam file
    # handle pvar file
    return bytes_written
end

function write_variant_record(io, xj::AbstractVector) # xj is the jth column of x
    N = length(xj)
    # track #4
    bytes_written = 0
    for xij in xj
        bytes_written += write(io, dosage_to_uint16(xij)) # 2 bytes per entry
    end
    return bytes_written
end

function dosage_to_uint16(xij::AbstractFloat, ploidy::Int=2)
    if isnan(xij)
        return int2bytes(65535, len=2)
    else
        return int2bytes(round(Int, xij/ploidy * 2^15), len=2)
    end
end
function dosage_to_uint16(::Missing, args...)
    return int2bytes(65535, len=2)
end

"""
int2bytes(x::Integer; len::Integer, little_endian::Bool)
    -> Vector{len, UInt8}

Convert an Integer `x` to a Vector{UInt8}
Options (not available for `x::BigInt`):
- `len` to define a minimum Vector lenght in bytes, result will show no leading
zero by default.
- set `little_endian` to `true` for a result in little endian byte order, result
in big endian order by default.
    julia> int2bytes(32974)
    2-element Array{UInt8,1}:
     0x80
     0xce
    julia> int2bytes(32974, len=4)
    4-element Array{UInt8,1}:
     0x00
     0x00
     0x80
     0xce
    julia> int2bytes(32974, little_endian=true)
    2-element Array{UInt8,1}:
     0xce
     0x80

# Source
https://github.com/roshii/BitConverter.jl/blob/master/src/BitConverter.jl
"""
function int2bytes(x::Integer; len::Integer=0, little_endian::Bool=true)
    result = reinterpret(UInt8, [hton(x)])
    i = findfirst(x -> x != 0x00, result)
    if len != 0
        i = length(result) - len + 1
    end
    result = result[i:end]
    if little_endian
        reverse!(result)
    end
    return result
end

"""
    bitstring2byte(s)

Parse a 8-digit bitstring (stored in Big endian) to an UInt8.

# Examples
+ bitstring2byte("00000001") = 0x01
+ bitstring2byte("01110001") = 0x71
"""
function bitstring2byte(s::AbstractString)
    @assert length(s) == 8
    return parse(UInt8, s, base=2)
end
