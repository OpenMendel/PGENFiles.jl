module PGEN
using BitIntegers
import Mmap: mmap
BitIntegers.@define_integers 24
const variant_type_lengths = Dict(
    0x00 => (4, 1), 0x01 => (4, 2), 0x02 => (4, 3), 0x03 => (4, 4),
    0x04 => (8, 1), 0x05 => (8, 2), 0x06 => (8, 3), 0x07 => (8, 4)
)
const bytes_to_UInt = Dict(1 => UInt8, 2 => UInt16, 3 => UInt24, 4 => UInt32, 8 => UInt64)
const mask_map = [0x01, 0x03, 0x00, 0x0f, 0x00, 0x00, 0x00, 0xff]
@inline ceil_int(x, y) = (x ÷ y) + (x % y != 0)
include("structs.jl")
include("header.jl")
end
