"""
    decode_single(input::AbstractVector{UInt8}; offset=zero(UInt))
Decode a single value of ULEB128 from `input`, starting with `offset + 1`-th byte.
"""
@inline function decode_single(input::Ptr{UInt8}; 
    offset::UInt=zero(UInt)
    )
    output = zero(UInt32)
    shift = 0x00
    for j in 1:5
        byte = unsafe_load(input, offset + j)#input[offset + j]
        output |= (UInt32(byte & 0x7f) << shift)
        if (byte & 0x80 == 0) 
            offset += j
            break
        end
        @assert j < 5 "incorrect encoding of a ULEB128 number"
        shift += 0x07
    end
    return output, offset
end

"""
    decode_multiple!(output::AbstractVector{UInt32}, input::AbstractVector{UInt8}; 
    offset=zero(UInt))

Decode `count` values of ULEB128 from `input`, starting with `offset + 1`-th byte.
"""
function decode_multiple!(output::Ptr{UInt32}, 
    input::Ptr{UInt8};
    count::Integer = length(output), 
    offset::UInt = zero(UInt)
    )
    for i in 1:count
        o = zero(UInt32)
        shift = 0x00
        for j in 1:5
            byte = unsafe_load(input, offset + j) #input[offset + j]
            o |= (UInt32(byte & 0x7f) << shift)
            if (byte & 0x80 == 0)
                offset += j 
                break
            end
            @assert j < 5 "incorrect encoding of a ULEB128 number"
            shift += 0x07
        end
        unsafe_store!(output, o, i) # output[i] = o
    end
    return output, offset
end

function size_n(input::AbstractVector{UInt8}, n::Integer, offset::UInt)
    r = zero(UInt32)
    for i in 1:n
        for j in 1:5
            byte = input[offset + j]
            if (byte & 0x80 == 0)
                r += j 
                offset += j
                break
            end
            @assert j < 5 "incorrect encoding of a ULEB128 number"
        end
    end
    r
end
