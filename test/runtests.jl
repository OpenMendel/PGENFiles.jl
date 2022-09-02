using PGENFiles
using Test, BGEN
const data = PGENFiles.datadir("bgen_example.16bits.pgen")
@testset "PGENFiles.jl" begin

function convert_gt(t::Type{T}, b::Bgen) where T <: Real
    n = BGEN.n_samples(b)
    p = BGEN.n_variants(b)

    # return arrays
    G = Matrix{t}(undef, n, p)
    Gchr = Vector{String}(undef, p)
    Gpos = Vector{Int}(undef, p)
    GsnpID = [String[] for _ in 1:p] # each variant can have >1 rsid, although we don't presently allow this
    Gref = Vector{String}(undef, p)
    Galt = [String[] for _ in 1:p] # each variant can have >1 alt allele, although we don't presently allow this

    # loop over each variant
    i = 1
    for v in BGEN.iterator(b; from_bgen_start=true)
        dose = first_allele_dosage!(b, v; T=t) # this reads REF allele as 1
        copyto!(@view(G[:, i]), dose)
        # store chr/pos/snpID/ref/alt info
        Gchr[i], Gpos[i] = chrom(v), pos(v)
        push!(GsnpID[i], rsid(v))
        ref_alt_alleles = alleles(v)
        length(ref_alt_alleles) > 2 && error("Marker $i of BGEN is not biallelic!")
        Gref[i] = ref_alt_alleles[1]
        push!(Galt[i], ref_alt_alleles[2])
        i += 1
        clear!(v)
    end

    return G, b.samples, Gchr, Gpos, GsnpID, Gref, Galt
end

function convert_gt(t::Type{T}, pfile::Pgen) where T <: Real
    n = PGENFiles.n_samples(pfile) |> Int
    p = PGENFiles.n_variants(pfile) |> Int

    # return arrays
    G = Matrix{t}(undef, n, p)

    # loop over each variant
    d = Vector{t}(undef, n)
    g = Vector{UInt8}(undef, n)
    g_ld = similar(g)
    for (j, v) in enumerate(PGENFiles.iterator(pfile))
        alt_allele_dosage!(d, g, pfile, v; genoldbuf=g_ld)
        v_rt = v.record_type & 0x07
        if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
            g_ld .= g
        end
        # store dosages
        G[:, j] .= d
    end

    return G
end

@testset "Header" begin
    p = PGENFiles.Pgen(data)
    h = p.header
    @test h.bits_per_variant_type == 8
    @test h.bytes_per_record_length == 2
    @test h.n_blocks == 1
    @test h.n_samples == 0x01f4
    @test h.n_variants == 0xc7
    @test h.provisional_reference == 0x01
    @test h.provisional_reference_flags === nothing
    @test h.allele_counts === nothing
    @test h.storage_mode == 0x10
    @test h.variant_block_offsets[1] == 0x0269
    @test h.variant_lengths[1] == 0x04a2
    @test h.variant_lengths[end] == 0x04a0
    @test h.variant_types[1] == 0x60
    @test h.variant_types[end-1] == 0x41
    @test h.bytes_per_sample_id == 0x02
end

@testset "Difflist" begin
    # Write an example difflist given in PGEN spec.
    # Note that an offset of 1 is added to sample indexes to make it 1-based.
    io = open("dummy", "w")

    write(io, 0x4f)
    
    write(io, 0x88)
    write(io, 0x13)
    write(io, 0x00)

    write(io, 0x88)
    write(io, 0xf5)
    write(io, 0x04)

    write(io, 0x3f)
    
    for i in 1:20
        write(io, 0x00)
    end
    
    for i in 1:77
        write(io, 0x88)
        write(io, 0x27)
    end
    close(io)

    d = read("dummy")
    dl, offset = PGENFiles.parse_difflist(d, UInt(0), 3, true)
    @test dl.len == 79
    @test all(dl.genotypes .== 0)
    @test dl.has_genotypes
    @test unsafe_load(dl.last_component_sizes, 1) == 0x3f
    @test dl.sample_id_bases[1] == 5000
    @test dl.sample_id_bases[2] == 325000
    #@test length(dl.sample_id_increments[]) == 154
    idx = Vector{UInt32}(undef, 64)
    idx_incr = Vector{UInt32}(undef, 64)
    PGENFiles.parse_difflist_sampleids!(idx, idx_incr, dl, 1)
    @test all(idx .== [5000 * i for i in 1:64] .+ 1)
    PGENFiles.parse_difflist_sampleids!(idx, idx_incr, dl, 2)
    @test all(idx[1:15] .== [5000 * (64 + i) for i in 1:15] .+ 1) # for idx 65..79
    @test all(idx[16:end] .== 0)

    rm("dummy", force=true)
end

@testset "dosage" begin
    # NOTE: First alleles in the BGEN file are encoded as alternate allele in the transformation.
    # Some record types are not covered by this test, LD-compressions in particular.
    # They have been tested on a private UK Biobank data file.
    b = Bgen(PGENFiles.datadir("example.16bits.bgen"))
    p = PGENFiles.Pgen(data)
    g_pgen = Array{UInt8}(undef, p.header.n_samples)
    g_pgen_ld = similar(g_pgen)
    d_pgen = Array{Float64}(undef, p.header.n_samples)
    for (v_bgen, v_pgen) in zip(BGEN.iterator(b), PGENFiles.iterator(p)) # 
        d_bgen = BGEN.first_allele_dosage!(b, v_bgen)
        PGENFiles.alt_allele_dosage!(d_pgen, g_pgen, p, v_pgen)      
        @test all(isapprox.(d_bgen, d_pgen; atol=5e-5, nans=true))
        PGENFiles.alt_allele_dosage!(d_pgen, g_pgen, p, v_pgen; genoldbuf=g_pgen_ld)  
        @test all(isapprox.(d_bgen, d_pgen; atol=5e-5, nans=true))
        v_rt = v_pgen.record_type & 0x07
        if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
            g_pgen_ld .= g_pgen
        end
    end
end

@testset "write PGEN" begin    
    # bitstring2byte function
    @test PGENFiles.bitstring2byte("01110001") == 0x71
    @test PGENFiles.bitstring2byte("11111111") == 0xff

    # bytes function (example at 2.2.6)
    @test PGENFiles.int2bytes(99325313, len=8) == [0x81, 0x95, 0xeb, 0x05, 0x00, 0x00, 0x00, 0x00]

    # dosage_to_uint16 function (example at end of 2.3.5)
    ploidy = 2
    @test PGENFiles.dosage_to_uint16(0.75, ploidy) == [0x00, 0x30]
    @test PGENFiles.dosage_to_uint16(1.5, ploidy) == [0x00, 0x60]

    # write_PGEN function
    bfile = Bgen(PGENFiles.datadir("example.16bits.bgen"))
    bgenG, nsamples, Gchr, Gpos, GsnpID, Gref, Galt = convert_gt(Float64, bfile)
    # pgenG = convert_gt(Float64, PGENFiles.Pgen(data))
    write_PGEN("test_pgen_write.pgen", bgenG)
    pgenG = convert_gt(Float64, PGENFiles.Pgen("test_pgen_write.pgen"))
    @test all(isapprox.(pgenG, bgenG; atol=5e-5, nans=true))
    rm("test_pgen_write.pgen", force=true)
end
end
