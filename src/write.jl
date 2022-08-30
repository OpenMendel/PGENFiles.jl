"""
    write_PGEN(pgen_filename, x, sampleID, variantID)

Saves numeric matrix `x` into a PGEN formatted file. We assume `x`
stores dosage genotype (i.e. x values are between 0 and 2). 
"""
function write_PGEN(
    pgen_filename::AbstractString, 
    x::AbstractMatrix,
    sampleID::AbstractVector,
    variantID::AbstractVector
    )
    n_samples, n_variants = size(x)
    # main PGEN file
    open(pgen_filename, "w") do io
        #
        # Construct header
        #
        # magic number
        write(io, 0x6c, 0x1b)
        # storage mode
        write(io, 0x10) # 2 bytes so far (note: 0-based indexing according to manual) 
        # data dimension
        write(io, UInt32(n_variants)) # 6 bytes
        write(io, UInt32(n_samples)) # 10 bytes
        # what is 11th byte?? 
        bits_per_record_type = 4
        bytes_per_record_length = 2
        write(io, 0x81) # following example in 2.2.6 for now
        # variant block offsets (i.e. start position for each block)
        n_blocks = PGENFiles.ceil_int(n_variants, 2^16) #Int(ceil(n_variants / 2 ^ 16))
        variant_offset = 12 + 8n_blocks
        variant_type_offset = Int(2^16 / (8 / bits_per_record_type))
        variant_length_offset = 2^16 * bytes_per_record_length
        variant_offset += (n_blocks - 1) * (variant_type_offset + variant_length_offset)
        last_block_variants = n_variants % 2^16
        variant_offset += Int(last_block_variants / (8 / bits_per_record_type)) + 
            last_block_variants * bytes_per_record_length # 99325313
        for b in n_blocks
            offset = b * variant_offset
            for x in bytes(variant_offset, len=8, little_endian=true)
                write(io, x)
            end
        end
        # variant record types
        #
        # Now handle variant records
        #
    end
    # handle pvar file
    # pvar_filename = joinpath(pgen_filename, )
    write_pvar(pvar_filneame)
end

function write_variant_record(io, xi::AbstractVector) # xi is a column of x
    N = length(xi)
    # track #3, assumes all samples have dosages
    bytes_written = 0
    for i in 1:div(N, 8)
        bytes_written += write(io, 0xff) # 0xff is UInt8 of 11111111
    end
    leftover = N % 8
    bytes_written += write(io, bitstring2byte("1"^leftover * "0"^(8 - leftover)))
    # track #4
    for xij in xi
        bytes_written += write(io, dosage_to_uint16(xij)) # 2 bytes per entry
    end
    return bytes_written
end

function dosage_to_uint16(xij::AbstractFloat, ploidy::Int=2)
    return bytes(round(Int, xij/ploidy * 2^15), len=2)
end

function write_pvar(pvar_filename)
    # todo
end

"""
    bytes(x::Integer; len::Integer, little_endian::Bool)
    -> Vector{len, UInt8}

Convert an Integer `x` to a Vector{UInt8}
Options (not available for `x::BigInt`):
- `len` to define a minimum Vector lenght in bytes, result will show no leading
zero by default.
- set `little_endian` to `true` for a result in little endian byte order, result
in big endian order by default.
    julia> bytes(32974)
    2-element Array{UInt8,1}:
     0x80
     0xce
    julia> bytes(32974, len=4)
    4-element Array{UInt8,1}:
     0x00
     0x00
     0x80
     0xce
    julia> bytes(32974, little_endian=true)
    2-element Array{UInt8,1}:
     0xce
     0x80

# Source
https://github.com/roshii/BitConverter.jl/blob/master/src/BitConverter.jl
"""
function bytes(x::Integer; len::Integer=0, little_endian::Bool=true)
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

Parse a 8-digit bitstring to an UInt8.

e.g. bitstring2byte("01110001") = 0x71
"""
function bitstring2byte(s::AbstractString)
    @assert length(s) == 8
    return parse(UInt8, s, base=2)
end
