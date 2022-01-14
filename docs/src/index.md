# PGEN.jl

Routines for reading compressed storage of genotyped or imputed markers

[*Genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) data with imputed markers are often saved in the [**PGEN format**](https://www.cog-genomics.org/plink/2.0/input#pgen) in `.pgen` file.
It can store both hard calls and imputed data, unphased genotypes and phased haplotypes. Each variant is compressed separately. This is the central data format for [PLINK 2](https://www.cog-genomics.org/plink/2.0/). 

## Format description



## Installation

This package requires Julia v1.6 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.


This package is registered in the default Julia package registry, and can be installed through standard package installation procedure: e.g., running the following code in Julia REPL.
```julia
using Pkg
pkg"add https://github.com/OpenMendel/PGEN.jl"
```



```julia
versioninfo()
```

    Julia Version 1.7.1
    Commit ac5cc99908 (2021-12-22 19:35 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin19.5.0)
      CPU: Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-12.0.1 (ORCJIT, skylake)


## Type `Pgen`

The type `Pgen` is the fundamental type for .pgen-formatted files. It can be created using the following line.




```julia
using PGEN
p = Pgen(PGEN.datadir("bgen_example.16bits.pgen")) ;
```

This example file is a PGEN file converted from a BGEN file.

Number of variants and samples in the file is accessible with the functions `n_variants()` and `n_samples()`. 


```julia
println(n_variants(p))
```

    199



```julia
println(n_samples(p))
```

    500


## `VariantIterator`

Genotype information of each variant is compressed separately in PGEN files. The offsets (starting points in pgen file) of each variant can be inferred from the header. A way to iterate over variants in a PGEN file is provided through a `VariantIterator` object, created by the function `iterator()`. 


```julia
v_iter = iterator(p);
```

One may check index, offset (starting point of each variant record), record type (how each variant record is comprssed and what type of information it has), and length of each variant using the following code:


```julia
for v in v_iter
    println("Variant $(v.index): offset $(v.offset), type 0x$(string(v.record_type, base=16))," * 
        " length $(v.length)")
end
```

    Variant 1: offset 617, type 0x60, length 1186
    Variant 2: offset 1803, type 0x60, length 1182
    Variant 3: offset 2985, type 0x41, length 1098
    Variant 4: offset 4083, type 0x60, length 1182
    Variant 5: offset 5265, type 0x60, length 1184
    Variant 6: offset 6449, type 0x60, length 1186
    Variant 7: offset 7635, type 0x44, length 1057
    Variant 8: offset 8692, type 0x61, length 1146
    Variant 9: offset 9838, type 0x61, length 1120
    Variant 10: offset 10958, type 0x60, length 1184
    Variant 11: offset 12142, type 0x61, length 1144
    Variant 12: offset 13286, type 0x40, length 1125
    Variant 13: offset 14411, type 0x41, length 1078
    Variant 14: offset 15489, type 0x40, length 1125
    Variant 15: offset 16614, type 0x41, length 1079
    Variant 16: offset 17693, type 0x66, length 1115
    Variant 17: offset 18808, type 0x60, length 1186
    Variant 18: offset 19994, type 0x44, length 1056
    Variant 19: offset 21050, type 0x60, length 1186
    Variant 20: offset 22236, type 0x46, length 1060
    Variant 21: offset 23296, type 0x61, length 1145
    Variant 22: offset 24441, type 0x61, length 1147
    Variant 23: offset 25588, type 0x44, length 1034
    Variant 24: offset 26622, type 0x41, length 1072
    Variant 25: offset 27694, type 0x41, length 1084
    Variant 26: offset 28778, type 0x40, length 1125
    Variant 27: offset 29903, type 0x61, length 1162
    Variant 28: offset 31065, type 0x44, length 1011
    Variant 29: offset 32076, type 0x41, length 1083
    Variant 30: offset 33159, type 0x40, length 1125
    Variant 31: offset 34284, type 0x61, length 1129
    Variant 32: offset 35413, type 0x44, length 1055
    Variant 33: offset 36468, type 0x61, length 1118
    Variant 34: offset 37586, type 0x60, length 1184
    Variant 35: offset 38770, type 0x40, length 1125
    Variant 36: offset 39895, type 0x60, length 1184
    Variant 37: offset 41079, type 0x61, length 1154
    Variant 38: offset 42233, type 0x40, length 1125
    Variant 39: offset 43358, type 0x60, length 1184
    Variant 40: offset 44542, type 0x41, length 1079
    Variant 41: offset 45621, type 0x41, length 1076
    Variant 42: offset 46697, type 0x61, length 1142
    Variant 43: offset 47839, type 0x60, length 1186
    Variant 44: offset 49025, type 0x61, length 1155
    Variant 45: offset 50180, type 0x41, length 1078
    Variant 46: offset 51258, type 0x40, length 1125
    Variant 47: offset 52383, type 0x41, length 1068
    Variant 48: offset 53451, type 0x41, length 1077
    Variant 49: offset 54528, type 0x40, length 1125
    Variant 50: offset 55653, type 0x60, length 1186
    Variant 51: offset 56839, type 0x40, length 1125
    Variant 52: offset 57964, type 0x41, length 1079
    Variant 53: offset 59043, type 0x41, length 1099
    Variant 54: offset 60142, type 0x60, length 1186
    Variant 55: offset 61328, type 0x40, length 1125
    Variant 56: offset 62453, type 0x40, length 1125
    Variant 57: offset 63578, type 0x40, length 1125
    Variant 58: offset 64703, type 0x40, length 1125
    Variant 59: offset 65828, type 0x40, length 1125
    Variant 60: offset 66953, type 0x61, length 1126
    Variant 61: offset 68079, type 0x60, length 1186
    Variant 62: offset 69265, type 0x40, length 1125
    Variant 63: offset 70390, type 0x44, length 1029
    Variant 64: offset 71419, type 0x40, length 1125
    Variant 65: offset 72544, type 0x60, length 1186
    Variant 66: offset 73730, type 0x60, length 1186
    Variant 67: offset 74916, type 0x40, length 1125
    Variant 68: offset 76041, type 0x41, length 1101
    Variant 69: offset 77142, type 0x61, length 1146
    Variant 70: offset 78288, type 0x60, length 1180
    Variant 71: offset 79468, type 0x41, length 1072
    Variant 72: offset 80540, type 0x44, length 1016
    Variant 73: offset 81556, type 0x60, length 1186
    Variant 74: offset 82742, type 0x44, length 1045
    Variant 75: offset 83787, type 0x40, length 1125
    Variant 76: offset 84912, type 0x40, length 1125
    Variant 77: offset 86037, type 0x41, length 1070
    Variant 78: offset 87107, type 0x60, length 1186
    Variant 79: offset 88293, type 0x61, length 1142
    Variant 80: offset 89435, type 0x44, length 1072
    Variant 81: offset 90507, type 0x41, length 1071
    Variant 82: offset 91578, type 0x60, length 1186
    Variant 83: offset 92764, type 0x41, length 1086
    Variant 84: offset 93850, type 0x40, length 1125
    Variant 85: offset 94975, type 0x61, length 1147
    Variant 86: offset 96122, type 0x41, length 1088
    Variant 87: offset 97210, type 0x60, length 1186
    Variant 88: offset 98396, type 0x40, length 1125
    Variant 89: offset 99521, type 0x44, length 1062
    Variant 90: offset 100583, type 0x40, length 1125
    Variant 91: offset 101708, type 0x41, length 1080
    Variant 92: offset 102788, type 0x46, length 1022
    Variant 93: offset 103810, type 0x40, length 1125
    Variant 94: offset 104935, type 0x60, length 1182
    Variant 95: offset 106117, type 0x40, length 1125
    Variant 96: offset 107242, type 0x40, length 1125
    Variant 97: offset 108367, type 0x61, length 1150
    Variant 98: offset 109517, type 0x41, length 1071
    Variant 99: offset 110588, type 0x60, length 1186
    Variant 100: offset 111774, type 0x60, length 1186
    Variant 101: offset 112960, type 0x40, length 1125
    Variant 102: offset 114085, type 0x40, length 1125
    Variant 103: offset 115210, type 0x61, length 1159
    Variant 104: offset 116369, type 0x60, length 1182
    Variant 105: offset 117551, type 0x60, length 1182
    Variant 106: offset 118733, type 0x60, length 1186
    Variant 107: offset 119919, type 0x46, length 1057
    Variant 108: offset 120976, type 0x41, length 1085
    Variant 109: offset 122061, type 0x61, length 1118
    Variant 110: offset 123179, type 0x60, length 1186
    Variant 111: offset 124365, type 0x61, length 1146
    Variant 112: offset 125511, type 0x40, length 1125
    Variant 113: offset 126636, type 0x41, length 1078
    Variant 114: offset 127714, type 0x40, length 1125
    Variant 115: offset 128839, type 0x41, length 1079
    Variant 116: offset 129918, type 0x44, length 1054
    Variant 117: offset 130972, type 0x60, length 1182
    Variant 118: offset 132154, type 0x46, length 1056
    Variant 119: offset 133210, type 0x60, length 1186
    Variant 120: offset 134396, type 0x44, length 1060
    Variant 121: offset 135456, type 0x61, length 1147
    Variant 122: offset 136603, type 0x61, length 1147
    Variant 123: offset 137750, type 0x46, length 1034
    Variant 124: offset 138784, type 0x41, length 1072
    Variant 125: offset 139856, type 0x41, length 1084
    Variant 126: offset 140940, type 0x40, length 1125
    Variant 127: offset 142065, type 0x61, length 1160
    Variant 128: offset 143225, type 0x46, length 1011
    Variant 129: offset 144236, type 0x41, length 1083
    Variant 130: offset 145319, type 0x60, length 1186
    Variant 131: offset 146505, type 0x41, length 1068
    Variant 132: offset 147573, type 0x46, length 1055
    Variant 133: offset 148628, type 0x61, length 1122
    Variant 134: offset 149750, type 0x60, length 1182
    Variant 135: offset 150932, type 0x40, length 1125
    Variant 136: offset 152057, type 0x60, length 1186
    Variant 137: offset 153243, type 0x61, length 1152
    Variant 138: offset 154395, type 0x40, length 1125
    Variant 139: offset 155520, type 0x60, length 1184
    Variant 140: offset 156704, type 0x61, length 1140
    Variant 141: offset 157844, type 0x41, length 1076
    Variant 142: offset 158920, type 0x61, length 1142
    Variant 143: offset 160062, type 0x60, length 1186
    Variant 144: offset 161248, type 0x61, length 1153
    Variant 145: offset 162401, type 0x41, length 1078
    Variant 146: offset 163479, type 0x40, length 1125
    Variant 147: offset 164604, type 0x41, length 1068
    Variant 148: offset 165672, type 0x41, length 1077
    Variant 149: offset 166749, type 0x40, length 1125
    Variant 150: offset 167874, type 0x60, length 1186
    Variant 151: offset 169060, type 0x40, length 1125
    Variant 152: offset 170185, type 0x41, length 1079
    Variant 153: offset 171264, type 0x61, length 1160
    Variant 154: offset 172424, type 0x60, length 1186
    Variant 155: offset 173610, type 0x40, length 1125
    Variant 156: offset 174735, type 0x40, length 1125
    Variant 157: offset 175860, type 0x40, length 1125
    Variant 158: offset 176985, type 0x40, length 1125
    Variant 159: offset 178110, type 0x60, length 1186
    Variant 160: offset 179296, type 0x61, length 1124
    Variant 161: offset 180420, type 0x60, length 1186
    Variant 162: offset 181606, type 0x40, length 1125
    Variant 163: offset 182731, type 0x46, length 1029
    Variant 164: offset 183760, type 0x40, length 1125
    Variant 165: offset 184885, type 0x40, length 1125
    Variant 166: offset 186010, type 0x60, length 1186
    Variant 167: offset 187196, type 0x60, length 1186
    Variant 168: offset 188382, type 0x41, length 1101
    Variant 169: offset 189483, type 0x61, length 1146
    Variant 170: offset 190629, type 0x60, length 1184
    Variant 171: offset 191813, type 0x41, length 1072
    Variant 172: offset 192885, type 0x46, length 1016
    Variant 173: offset 193901, type 0x60, length 1184
    Variant 174: offset 195085, type 0x46, length 1045
    Variant 175: offset 196130, type 0x60, length 1186
    Variant 176: offset 197316, type 0x40, length 1125
    Variant 177: offset 198441, type 0x41, length 1070
    Variant 178: offset 199511, type 0x60, length 1184
    Variant 179: offset 200695, type 0x61, length 1142
    Variant 180: offset 201837, type 0x46, length 1072
    Variant 181: offset 202909, type 0x61, length 1132
    Variant 182: offset 204041, type 0x60, length 1186
    Variant 183: offset 205227, type 0x41, length 1086
    Variant 184: offset 206313, type 0x40, length 1125
    Variant 185: offset 207438, type 0x41, length 1088
    Variant 186: offset 208526, type 0x41, length 1088
    Variant 187: offset 209614, type 0x60, length 1184
    Variant 188: offset 210798, type 0x40, length 1125
    Variant 189: offset 211923, type 0x46, length 1062
    Variant 190: offset 212985, type 0x60, length 1186
    Variant 191: offset 214171, type 0x41, length 1080
    Variant 192: offset 215251, type 0x44, length 1022
    Variant 193: offset 216273, type 0x40, length 1125
    Variant 194: offset 217398, type 0x60, length 1186
    Variant 195: offset 218584, type 0x40, length 1125
    Variant 196: offset 219709, type 0x40, length 1125
    Variant 197: offset 220834, type 0x61, length 1152
    Variant 198: offset 221986, type 0x41, length 1071
    Variant 199: offset 223057, type 0x60, length 1184


More information on each variant is available in the attached `.pvar` file. 

## Genotypes and dosages

Genotypes of each variant is available through the function `get_genotypes()` or `get_genotypes!()`. For example, to obtain the genotypes of the first variant, one may do:


```julia
v = first(v_iter)
g, data, offset = get_genotypes(p, v)
```




    (UInt8[0x03, 0x00, 0x00, 0x01, 0x00, 0x03, 0x01, 0x00, 0x03, 0x03  …  0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x03, 0x01], UInt8[0x43, 0x1c, 0xff, 0x14, 0xc7, 0x0f, 0x00, 0x30, 0x01, 0x04  …  0xdc, 0x03, 0xd3, 0x42, 0x9e, 0x03, 0x07, 0x79, 0x4b, 0x3f], 0x000000000000007d)



`g` stores the genotypes, `data` is the variant record for `v`, and `offset` indicates where the track for genotypes ended. Encoding for `g` is as following:

| genotype code | genotype category | 
|:---:|:---:|
| `0x00` | homozygous REF | 
| `0x01` | heterozygous REF-ALT |
| `0x02` | homozygous ALT |
| `0x03` | missing |

To avoid array allocations for iterative genotype extraction, one may preallocate `g` and reuse it. As some of the variants are LD-compressed, an additional genotype buffer to keep the genotypes for the most recent non-LD-compressed variant may be desirable (`g_ld`). If `g_ld` is not provided, it will parse the genotypes of the most recent non-LD-compressed variant (stored in an internal dictionary) first.

For example:


```julia
g = Vector{UInt8}(undef, n_samples(p))
g_ld = similar(g)
for v in v_iter
    get_genotypes!(g, p, v; ldbuf=g_ld)
    v_rt = v.record_type & 0x07
    if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
        g_ld .= g
    end
    
    # do someting with genotypes in `g`...
end
```

Similarly, ALT allele dosages are available through the function `alt_allele_dosage()` and `alt_allele_dosage!()`. As genotype information is required to retrieve dosages, space for genotypes are also required for `alt_allele_dosage!()`. These functions return dosages, parsed genotypes, and endpoint of the dosage information in the current variant record.

To obtain the dosages of the first variant: 


```julia
v = first(v_iter)
d, g, offset = alt_allele_dosage(p, v)
```




    (Float32[NaN, 0.06427002, 0.08441162, 0.98254395, 0.08843994, 0.14111328, 1.0733032, 0.054138184, 0.10858154, 0.12310791  …  0.029785156, 0.9661255, 0.00079345703, 1.0126343, 0.042663574, 0.060302734, 1.0441284, 0.056518555, 1.8910522, 0.98895264], UInt8[0x03, 0x00, 0x00, 0x01, 0x00, 0x03, 0x01, 0x00, 0x03, 0x03  …  0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x03, 0x01], 0x00000000000004a2)



Missing value is represented by a `NaN`. Code for a typical GWAS application should look like:


```julia
d = Vector{Float32}(undef, n_samples(p))
g = Vector{UInt8}(undef, n_samples(p))
g_ld = similar(g)
for v in v_iter
    alt_allele_dosage!(d, g, p, v; genoldbuf=g_ld)
    v_rt = v.record_type & 0x07
    if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
        g_ld .= g
    end
    
    # do someting with dosage values in `d`...
end
```

## Speed

The current PGEN package can read in ~2000 variants / second for UK Biobank data, which is about 4x faster than reading in BGEN-formatted data. 
