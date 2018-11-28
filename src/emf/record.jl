
# EMF records constitue a single column in the alignment



"""
Metadata for one row in the multiple sequence alignment.
"""
mutable struct RowInterval
    data::Vector{UInt8}
    organism::UnitRange{Int}
    chrom::UnitRange{Int}
    first::UnitRange{Int}
    last::UnitRange{Int}
    strand::UnitRange{Int}
    chrlen::UnitRange{Int}

    function RowInterval()
        return new(UInt8[], 1:0, 1:0, 1:0, 1:0, 1:0, 1:0)
    end
end


"""
A record represents one column in the multiple sequence alignment.
"""
mutable struct Record
    score_type_data::Vector{UInt8}
    score_type::UnitRange{Int}

    tree_data::Vector{UInt8}
    tree::UnitRange{Int}

    rows::Vector{RowInterval}
    nrows::Int

    column_data::Vector{UInt8}
    column::UnitRange{Int}
    score::UnitRange{Int}

    function Record()
        return new(
            UInt8[], 1:0,
            UInt8[], 1:0,
            RowInterval[], 0,
            UInt8[], 1:0, 1:0)
    end
end


function nrows(rec::Record)
    return rec.nrows
end


function row(rec::Record, i::Integer)
    @assert 1 <= i <= rec.nrows
    return rec.rows[i]
end


function organism(row::RowInterval)
    return String(row.data[row.organism])
end


function chrom(row::RowInterval)
    return String(row.data[row.chrom])
end


function Base.first(row::RowInterval)
    return BioCore.RecordHelper.unsafe_parse_decimal(Int, row.data, row.first)
end


function Base.last(row::RowInterval)
    return BioCore.RecordHelper.unsafe_parse_decimal(Int, row.data, row.last)
end


function strand(row::RowInterval)
    if length(row.strand) == 2 &&
            row.data[first(row.strand)] == UInt8('-') &&
            row.data[last(row.strand)] == UInt8('1')
        return -1
    elseif length(row.strand) == 1 && row.data[first(row.strand)] == UInt8('1')
        return 1
    else
        error("Cannot parse strand \"", String(row.data[row.strand]), "\"")
    end
end


function chrlen(row::RowInterval)
    return BioCore.RecordHelper.unsafe_parse_decimal(Int, row.data, row.chrlen)
end


function score_type(rec::Record)
    return String(rec.score_type_data[rec.score_type])
end


function tree(rec::Record)
    return String(rec.tree_data[rec.tree])
end


function score(rec::Record)
    if length(rec.score) == 1 && rec.column_data[first(rec.score)] == UInt8('.')
        return nothing
    else
        hasvalue, val = ccall(:jl_try_substrtod,
            Tuple{Bool, Float64}, (Ptr{UInt8},Csize_t,Csize_t),
            rec.column_data, first(rec.score)-1, length(rec.score))
        if !hasvalue
            error("Invalid column score in EMF alignment column \"",
                String(rec.column_data[rec.score]), "\"")
        end
        return val
    end
end


function column(rec::Record)
    for i in rec.column
        if rec.column_data[i] == UInt8('.')
            rec.column_data[i] = 'N'
        end
    end
    return DNASequence(rec.column_data, first(rec.column), last(rec.column))
end


function initialize!(record::Record)
    record.score_type = 1:0
    record.tree = 1:0
    record.nrows = 0
    record.column = 1:0
    record.score = 1:0

    return record
end

function BioCore.isfilled(record::Record)
    return !isempty(record.column)
end
