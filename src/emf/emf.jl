# EMF (Ensembl Multi Format)
# ==========================

module EMF

using Automa
using Automa.RegExp: @re_str
using BioCore
using GenomicFeatures: Interval
using BioSequences
using BufferedStreams

include("record.jl")
include("reader.jl")

end
