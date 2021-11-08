module PGEN
using BitIntegers
import Mmap: mmap
BitIntegers.@define_integers 24
const bytes_to_UInt = Dict(1 => UInt8, 2 => UInt16, 3 => UInt24, 4 => UInt32, 8 => UInt64)
const mask_map = [0x01, 0x03, 0x00, 0x0f, 0x00, 0x00, 0x00, 0xff]
include("structs.jl")
include("header.jl")
end
