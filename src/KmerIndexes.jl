module KmerIndexes

export KmerIndex, count, locate

using Bio.Seq

"""
Dense k-mer index.
"""
immutable KmerIndex{K,T}
    data::Vector{Vector{T}}
end

function KmerIndex{K}(::Type{DNAKmer{K}}, seq::DNASequence, step::Integer=1)
    data = Vector{Vector{Int}}(4^K)
    for idx in 1:endof(data)
        data[idx] = Vector{Int}()
    end
    for (pos, kmer) in each(DNAKmer{K}, seq, step)
        idx = UInt64(kmer)
        push!(data[idx+1], pos)
    end
    return KmerIndex{K,Int}(data)
end

function Base.count{K}(kmer::DNAKmer{K}, index::KmerIndex{K})
    return length(locate(kmer, index))
end

function locate{K}(kmer::DNAKmer{K}, index::KmerIndex{K})
    idx = UInt64(kmer)
    return index.data[idx+1]
end

end # module
