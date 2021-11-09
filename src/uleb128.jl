@inline function decode_single(input::AbstractVector{UInt8}; 
    offset::Base.RefValue{UInt}=Ref(zero(UInt))
    )
    output = zero(UInt32)
    shift = 0x00
    offset_ = offset[]
    for j in 1:5
        byte = input[offset_ + j]
        output |= (UInt32(byte & 0x7f) << shift)
        if (byte & 0x80 == 0) 
            offset_ += j
            break
        end
        if j == 5
            @error "incorrect encoding of a ULEB128 number"
        end
        shift += 0x07
    end
    offset[] = offset_
    return output
end

function decode_multiple!(output::AbstractVector{UInt32}, 
    input::AbstractVector{UInt8};
    count::Integer = length(output), 
    offset::Base.RefValue{UInt} = Ref(zero(UInt))
    )
    offset_ = offset[]
    for i in 1:count
        o = zero(UInt32)
        shift = 0x00
        for j in 1:5
            byte = input[offset_ + j]
            o |= (UInt32(byte & 0x7f) << shift)
            if (byte & 0x80 == 0)
                offset_ += j 
                break
            end
            if j == 5
                @error "incorrect encoding of a ULEB128 number"
            end
            shift += 0x07
        end
        offset[] = offset_
        output[i] = o
    end
    return output
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
            if j == 5
                @error "incorrect encoding of a ULEB128 number"
            end
        end
    end
    r
end
