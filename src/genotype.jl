function get_genotypes!(buf::Vector{UInt8}, p::Pgen, v::Variant)

    buf
end

function get_genotypes(p::Pgen, varidx::Int)
    buf = Vector{UInt8}(undef, p.header.n_samples)
    get_genotypes!(buf, p, varidx)
end