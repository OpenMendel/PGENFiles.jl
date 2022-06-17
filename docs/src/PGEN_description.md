# The PGEN format

Content on this page is based on the [draft specification](https://github.com/chrchang/plink-ng/raw/master/pgen_spec/pgen_spec.pdf), distributed under GPLv3. 

## Introduction

The PGEN format is the central file format for genomic data in PLINK 2. 

- PLINK 1’s binary genotype file format (the BED format, can be read using [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl))
    - Simple, compact, and supports direct computation on the packed data representation. Thanks to these properties, it continues to be widely used more than a decade after it was designed.
    - Limitation: can only represent unphased biallelic hard-called genotypes.
        - suboptimal for GWASes which tend to benefit from inclusion of imputed dosages and more sophicated handling of multiallelic variants
        - cannot represent phase information for workflows like investigation of compound heterozygosity, imputation-related data management, etc.

The widely-used binary genotype formats which addresses limitations above include BCF format and BGEN format, but they do not support direct computation on packed data, impossible to match efficiency of PLINK 1.9. 
    
Hence, PLINK 2 decided to introduce a new binary genotype file format, the PGEN format. 

- Backward-compatible to BED format
- Can represent phased, multiallelic, and dosage data in a manner that better support ["compressive genomics"](https://www.nature.com/articles/nbt.2241)
- Incorporates "[SNPack](https://sysbiobig.dei.unipd.it/software/?q=software#SNPack)-style" genotype compression, reducing file sizes by 80+% with negligible encoding and decoding cost (and supporting some direct computation on the compressed representation)

Not as simple as PLINK 1 format, but now it includes open-source internal library (pgenlib) to read and write the format. 

Also introduced are: 
- PSAM format, an extension of `.fam` format
    - Stores categorical and other phenotype/covariates
- PVAR format, an extension of `.bim` format.
    - Stores all header and variant-specific information. 
    - Designed so that "sites-only VCF" files are directly valid PVAR files.

## PGEN format


Binary format capable of representing mixed-phase, multiallelic, mixed-hardcall/dosage/missing genotype data.

- A PLINK 1 variant-major .bed file is grandfathered in as a valid PGEN file. Simple to handle with the PGEN format definition.
- PGEN(+PVAR) is designed to interoperate with, not replace, VCF/BCF. 
    - PGEN cannot represent: read depths, quality scores, or biallelic genotype probability triplets, or triploid genotypes. 
    - It specializes on the subset of the VCF format which is relevant to PLINK’s function. 
    - Fast VCF ↔ PGEN conversion in PLINK 2.

### File organization

- Header: information to enable random access to the variant record: e.g., record types and record length of each variant.
    - Here, record type means how the genotype is compressed, if it contains phase and dosage information, etc. 
- A sequence of variant records.


A variant record’s main data track can be “LD-compressed” (LD = linkage disequilibrium):
- Most recent non-LD-compressed variant record and only storing genotype-category differences from it. 
    - The only type of inter-record dependency, 
- Record type and size information in the header, and the genotypes from the latest non-LD-compressed variant is enough to decode genotypes of each variant sequentially.

<!--
Three fixed-width storage modes are defined (covering basic unphased biallelic genotypes, unphased dosages, and phased dosages) which don’t have this limitation, and are especially straightforward
to read and write; but they don’t benefit from PGEN’s low-overhead genotype compression. A future version of this
specification may add a way to store most header information in a separate file, so that sequential reading, sequential
writing, and genotype compression are simultaneously possible (at the cost of more annoying file management).-->

### Header


- Magic number: `0x6c 0x1b`. 
- Storage mode
    - `0x01`: PLINK 1 BED format. Supported in `SnpArrays.jl`.
    - `0x02`: the simplest PLINK 2 fixed-width format for unphased genotypes. Difference from `0x01` are header and genotype encoding rule.
    - `0x03`: fixed-width unphased dosage
    - `0x04`: fixed-width phased dosage
    - __`0x10`__: standard variable-width format. Vast majority of the PLINK 2 files will be in this mode. __Currently, only this mode is supported in `PGEN.jl`__. 
- Dataset dimensions, header body formatting
    - number of variants, samples, bits per record type, bytes per allele counts (for multiallelic variants), if reference alleles are provisional, etc.
- Variant block offsets
    - Where each block of $2^{16}$ = 65,536 variant records begin. i.e. starting point of variant 1, 65_537, 131_073, ... 
- Main header body
    - Packed array of $2^{16}$ record types, lengths, etc. 
    
e.g., random access to 65540-th (65536 + 4) variant can be achieved by scannig for the 2nd entry of Variant block offsets and then scanning the first four entries of main header body. The starting point of the variant is calculated by the start of the second variant block plus first three variant record lengths. 

### Variant record

Each variant record starts with the main track for unphased biallelic hard-call genotypes, followed by the ten optional tracks:
1. Multiallelic hard-calls
2. Hardcall-phase information
3. Biallelic dosage existence
4. Biallelic dosage values
5. Multiallelic dosage existence
6. Multiallelic dosage values
7. Biallelic phased-dosage existence
8. Biallelic phased-dosage values
9. Multiallelic phased-dosage existence
10. Multiallelic phased-dosage values


### Difflists
Many genotypes and dosages are compressed in a __difflist__. It is designed to represent a sparse list of differences from something else. It does so in a manner that is compact, and supports fast checking of whether a specific sample ID is in the list. Struct for difflist is in the struct `DiffList`. 

### Main track

Each genotype is represented in two-bit little-endian ordering: e.g., the for the two bytes of `0x1b 0xd8` for 8 samples:

```
byte 1         byte 2
0x1b           0xd8
00 01 10 11    11 01 10 00
s4 s3 s2 s1    s8 s7 s6 s5
```

| Sample index (1-based) | genotype category | 
|:---:|:---:|
| 1 | `0b11` |
| 2 | `0b10` |
| 3 | `0b01` |
| 4 | `0b00` |
| 5 | `0b00` |
| 6 | `0b10` |
| 7 | `0b01` |
| 8 | `0b11` |

| genotype category | PLINK 1 | PLINK 2 | 
|:---:|:---:|:---:|
| 0 = `0b00` = `0x00` | homozygous A1 | homozygous REF | 
| 1 = `0b01` = `0x01` | missing | heterozygous REF-ALT |
| 2 = `0b10` = `0x02` | heterozygouus A1-A2 | homozygous ALT |
| 3 = `0b11` = `0x03` | homozygous A2 | missing |
 
 - A1: First allele listed in PLINK 1 bim file
 - A2: Second allele listed in PLINK 1 bim file
 - REF: Reference allele
 - ALT: Alternate allele
 
In PLINK 1, A1 was often ALT, and A2 was often REF. However, this was not set in stone. In UK Biobank data, A1 is REF and A2 is ALT. 

Seven record types are supported, represented by the bottom three bits of record type:

- `0`: no compression.
- `1`: “1-bit” representation. This starts with a byte indicating what the two most common categories are (value 1: categories 0 and 1; 2: 0 and 2; 3: 0 and 3; 5: 1 and 2; 6: 1 and 3; 9: 2 and 3); followed by a bitarray describing which samples are in the higher-numbered category; followed by a difflist with all (sample ID, genotype category value) pairs for the two less common categories.
- `2`: LD-compressed. A difflist with all (sample ID, genotype category value) pairs for samples in different categories than they were in in the previous non-LD-compressed variant. The first variant of a variant block (i.e. its index is congruent to 0 mod $2^{16}$) cannot be LD-compressed.
- `3`: LD-compressed, inverted. A difflist with all (sample ID, inverted genotype value) pairs for samples in different categories than they would be in the previous non-LD-compressed variant after inversion (categories 0 and 2 swapped). I.e. decoding can be handled in the same way as for variant record type 2, except for a final inversion step applied after the difflist contents have been patched in. This addresses spots where the reference genome is “wrong” for the population of interest.
- `4`: Difflist with all (sample ID, genotype category value) pairs for samples outside category 0.
- ~~`5`: Reserved for future use. (When all samples appear to be in category 1, that usually implies a systematic variant calling problem.)~~
- `6`: Difflist with all (sample ID, genotype category value) pairs for samples outside category 2.
- `7`: Difflist with all (sample ID, genotype category value) pairs for samples outside category 3

### Multiallelic hardcalls
Exists if the 4th bit of variant record type is set. Based on the main track, it defines a "patch set" in the form of difflist which sample has alternate allele other than "ALT1". 

###  Phased heterozygous hard-calls
Exists if the 5th bit of variant record type is set. Stores whether each heterozygous call is phased, and if phased, what the phase is. "`0|1`" or "`1|0`". PGEN does not distinguish "`0|0`" from "`0/0`", and "`1|1`" from "`1/1`". 

### Dosages
Dosages are stored in 16-bit integers (`UInt16`). `0x0000`...`0x8000`($2^{15}$) represent diploid ALT allele dosage values between `0.0`..`2.0`. `0xffff` represents missing value. Three record types are supported, based on 6th and 7th bits of record type. Dosages are required to be consistent with hard-calls (should be close enough from genotype).

- 6th bit is set and 7th bit is clear: Track 3 (Biallelic dosage existence) is a difflist indicating which samples have dosage information. 
- 6th bit is clear and 7th bit is set: Track 3 is omitted and Track 4 (Biallelic dosage values) has an entry for every single sample.
- 6th bit and 7th bit are both set: Track 3 is a BitArray indicating dosage for which sample exists. 

Samples without dosage values are assumed to have dosage level identical to their respective genotypes.
