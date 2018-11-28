# EMF Reader
# ==========

# TODO: parameterize by alphabet type???
mutable struct Reader <: BioCore.IO.AbstractReader
    state::BioCore.Ragel.State

    # file-wide metadata

    # "compara", "resequencing", or "gene_alignment". Currently only "compara",
    # i.e. multi-way whole genome alignments are supported.
    format::Union{Nothing, String}

    date::Union{Nothing, String}
    release::Union{Nothing, String}

    function Reader(input::BufferedInputStream)
        return new(
            BioCore.Ragel.State(body_machine.start_state, input),
            nothing, nothing, nothing)
    end
end



"""
    EMF.Reader(input::IO)

Create a data reader of the EMF file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO)
    return Reader(BufferedInputStream(input))
end

function BioCore.IO.stream(reader::Reader)
    return reader.state.stream
end


function Base.eltype(::Type{Reader})
    return Record
end


const body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt
    any = Automa.RegExp.any

    hspace = re"[ \t\v]+"

    directive = re"##[^\r\n]*"
    # directive.actions[:enter] = [:anchor]
    # directive.actions[:exit] = [:directive]

    comment = re"#([^#\r\n][^\r\n]*)?"

    seq_organism = re"[^ \t\v]*"
    seq_organism.actions[:enter] = [:mark]
    seq_organism.actions[:exit] = [:seq_organism]

    seq_chr = re"[^ \t\v]*"
    seq_chr.actions[:enter] = [:mark]
    seq_chr.actions[:exit] = [:seq_chr]

    seq_first = re"[0-9]+"
    seq_first.actions[:enter] = [:mark]
    seq_first.actions[:exit] = [:seq_first]

    seq_last = re"[0-9]+"
    seq_last.actions[:enter] = [:mark]
    seq_last.actions[:exit] = [:seq_last]

    seq_strand = re"\-?1"
    seq_strand.actions[:enter] = [:mark]
    seq_strand.actions[:exit] = [:seq_strand]

    seq_chrlen = re"[0-9]+"
    seq_chrlen.actions[:enter] = [:mark]
    seq_chrlen.actions[:exit] = [:seq_chrlen]

    seq = cat(
        "SEQ",
        hspace, seq_organism,
        hspace, seq_chr,
        hspace, seq_first,
        hspace, seq_last,
        hspace, seq_strand,
        opt(cat(hspace, "(chr_length=", seq_chrlen, ")")))
    seq.actions[:enter] = [:seq_enter]
    seq.actions[:exit] = [:seq_exit]

    score_type_data = re"[^\r\n]*"
    score_type_data.actions[:enter] = [:anchor]
    score_type_data.actions[:exit] = [:score_type_data]
    score_type = cat("SCORE", hspace, score_type_data)

    tree_data = re"[^\r\n]*"
    tree_data.actions[:enter] = [:anchor]
    tree_data.actions[:exit] = [:tree_data]
    tree = cat("TREE", hspace, tree_data)

    id_data = re"[^\r\n]*"
    # id_data.actions[:enter] = [:anchor]
    # id_data.actions[:exit] = [:id_data]
    id = cat("ID", hspace, id_data)

    # ensembl uses '.' to represent an unknown nucleotide
    letters = re"[A-Za-z*\-\.]+"
    letters.actions[:enter] = [:mark]
    letters.actions[:exit]  = [:letters]

    score = re"[ -~]*[0-9][ -~]*|\."
    score.actions[:enter] = [:mark]
    score.actions[:exit] = [:score]

    column = cat(letters, hspace, score)
    column.actions[:enter] = [:anchor]
    column.actions[:exit] = [:column]

    blank = re"[ \t]*"

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]
        cat(opt('\r'), lf)
    end

    block = cat(
        rep(cat(alt(seq, score_type, tree, id, directive, comment, blank), newline)),
        "DATA", newline,
        rep(cat(column, newline)),
        "//", newline)
    block.actions[:enter] = [:block_enter]

    body = rep(block)

    # return map(Automa.compile, (body,))
    return Automa.compile(body)
end)()


const actions = Dict(
    # :directive => quote
    #     # If we want to read particular directives, they should be added
    #     # to the grammar
    # end,

    :seq_organism => :(record.rows[record.nrows+1].organism = (mark:p-1) .- offset),
    :seq_chr      => :(record.rows[record.nrows+1].chrom    = (mark:p-1) .- offset),
    :seq_first    => :(record.rows[record.nrows+1].first    = (mark:p-1) .- offset),
    :seq_last     => :(record.rows[record.nrows+1].last     = (mark:p-1) .- offset),
    :seq_strand   => :(record.rows[record.nrows+1].strand   = (mark:p-1) .- offset),
    :seq_chrlen   => :(record.rows[record.nrows+1].chrlen   = (mark:p-1) .- offset),

    :seq_enter => quote
        while record.nrows+1 > length(record.rows)
            push!(record.rows, RowInterval())
        end
        BioCore.ReaderHelper.anchor!(stream, p); offset = p - 1
    end,

    :seq_exit => quote
        BioCore.ReaderHelper.resize_and_copy!(
            record.rows[record.nrows+1].data, data,
            BioCore.ReaderHelper.upanchor!(stream):p-1)
        record.nrows += 1
    end,

    :tree_data => quote
        BioCore.ReaderHelper.resize_and_copy!(
            record.tree_data, data,
            BioCore.ReaderHelper.upanchor!(stream):p-1)
        record.tree = (offset+1:p-1) .- offset
    end,

    :score_type_data => quote
        BioCore.ReaderHelper.resize_and_copy!(
            record.score_type_data, data,
            BioCore.ReaderHelper.upanchor!(stream):p-1)
        record.score_type = (offset+1:p-1) .- offset
    end,

    :letters => :(record.column = (mark:p-1) .- offset),
    :score   => :(record.score  = (mark:p-1) .- offset),

    :column    => quote
        BioCore.ReaderHelper.resize_and_copy!(
            record.column_data, data,
            BioCore.ReaderHelper.upanchor!(stream):p-1)
        found_record = true
        @escape
    end,

    :block_enter => quote
        for row in record.rows
            row.organism = 1:0
            row.chrom = 1:0
            row.first = 1:0
            row.last = 1:0
            row.strand = 1:0
            row.chrlen = 1:0
        end
    end,

    :countline => :(linenum += 1),
    :mark      => :(mark = p),
    :anchor    => :(BioCore.ReaderHelper.anchor!(stream, p); offset = p - 1)
)


eval(
    BioCore.ReaderHelper.generate_read_function(
        Reader,
        body_machine,
        quote
            mark = offset = 0
        end,
        actions
    )
)
