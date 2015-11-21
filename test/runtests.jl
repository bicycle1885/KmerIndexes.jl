using KmerIndexes
using BaseTestNext
using Bio.Seq

@testset "short sequence" begin
    seq = dna"ACGTAGACGT"
    index = KmerIndex(DNAKmer{4}, seq)
    @test count(kmer(dna"ACGT"), index) === 2
    @test locate(kmer(dna"ACGT"), index) == [1, 7]
    @test count(kmer(dna"AGAC"), index) === 1
    @test locate(kmer(dna"AGAC"), index) == [5]
end
