using PGEN
using Test, BGEN
const data = PGEN.datadir("bgen_example.16bits.pgen")
@testset "PGEN.jl" begin

@testset "Header" begin
    p = PGEN.Pgen(data)
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
    dl, offset = PGEN.parse_difflist(d, UInt(0), 3, true)
    @test dl.len == 79
    @test all(dl.genotypes .== 0)
    @test dl.has_genotypes
    @test dl.last_component_sizes[][1] == 0x3f
    @test dl.sample_id_bases[][1] == 5000
    @test dl.sample_id_bases[][2] == 325000
    @test length(dl.sample_id_increments[]) == 154
    idx = Vector{UInt32}(undef, 64)
    idx_incr = Vector{UInt32}(undef, 64)
    PGEN.parse_difflist_sampleids!(idx, idx_incr, dl, 1)
    @test all(idx .== [5000 * i for i in 1:64] .+ 1)
    PGEN.parse_difflist_sampleids!(idx, idx_incr, dl, 2)
    @test all(idx[1:15] .== [5000 * (64 + i) for i in 1:15] .+ 1) # for idx 65..79

    rm("dummy", force=true)
end

@testset "dosage" begin
    # NOTE: First alleles in the BGEN file are encoded as alternate allele in the transformation.
    # Some record types are not covered by this test, LD-compressions in particular.
    # They have been tested on a private UK Biobank data file.
    b = Bgen(PGEN.datadir("example.16bits.bgen"))
    p = PGEN.Pgen(data)
    g_pgen = Array{UInt8}(undef, p.header.n_samples)
    g_pgen_ld = similar(g_pgen)
    d_pgen = Array{Float64}(undef, p.header.n_samples)
    for (v_bgen, v_pgen) in zip(BGEN.iterator(b), PGEN.iterator(p)) # 
        d_bgen = BGEN.ref_allele_dosage!(b, v_bgen)
        PGEN.alt_allele_dosage!(d_pgen, g_pgen, p, v_pgen)      
        @test all(isapprox.(d_bgen, d_pgen; atol=5e-5, nans=true))
        PGEN.alt_allele_dosage!(d_pgen, g_pgen, p, v_pgen; genoldbuf=g_pgen_ld)  
        @test all(isapprox.(d_bgen, d_pgen; atol=5e-5, nans=true))
        v_rt = v_pgen.record_type & 0x07
        if v_rt != 0x02 && v_rt != 0x03 # non-LD-compressed. See Format description.
            g_pgen_ld .= g_pgen
        end
    end
end
end
