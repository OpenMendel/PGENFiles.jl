@inline function decode_single(input::AbstractVector{UInt8}; offset::UInt=zero(UInt))
    output = zero(UInt32)
    shift = 0x00
    for j in 1:5
        byte = input[offset + j]
        output |= (UInt32(byte & 0x7f) << shift)
        if (byte & 0x80 == 0) 
            break
        end
        shift += 0x07
    end
    return output
end

function decode_multiple!(output::AbstractVector{UInt32}, 
    input::AbstractVector{UInt8};
    count::UInt8 = length(output), 
    offset::UInt = zero(UInt))
    for i in 1:count
        o = zero(UInt32)
        shift = 0x00
        for j in 1:5
            byte = input[offset + j]
            o |= (UInt32(byte & 0x7f) << shift)
            if (byte & 0x80 == 0)
                offset += j 
                break
            end
            shift += 0x07
        end
        output[i] = o
    end
    return output
end