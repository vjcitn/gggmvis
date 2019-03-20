#' demonstrate ggplot-based gene model concept
#' @importFrom GenomicFeatures exons
#' @param sym character(1) gene symbol
#' @param resource an instance of \code{\link{EnsDb-class}}
#' @param yval numeric(1) vertical position of segments
#' @param arrmm numeric(1) number of millimeters for arrowhead segments
#' @param viewtype character(1) either "exons" or "transcripts"
#' @param \dots passed to GenomicFeatures::exons
#' @return ggplot instance
#' @note It is known that arrow heads do not propagate to ggplotly
#' displays as of Aug 2018.
#' @examples
#' requireNamespace("EnsDb.Hsapiens.v75")
#' edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
#' ggvisForSymbol("ORMDL3", arrmm=1.3, resource=edb) + ylab(" ") + 
#'       theme(axis.text.y = element_blank()) + ylim(-1,3)
#' @export
ggvisForSymbol = function (sym, resource = EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79, 
    columnsKept = c("gene_id", "tx_id"), yval = 1, arrmm=1.5, viewtype="transcripts", ...) 
{
    exs = GenomicFeatures::exons(resource, filter = SymbolFilter(sym), columns = columnsKept, 
        ...)
    if (viewtype == "exons") exs = unique(exs)
    rd = reduce(exs)
    fo = findOverlaps(rd, exs)
    gr = split(subjectHits(fo), queryHits(fo))
    pp = function(n) (seq_len(n)-1)/n
    st = start(exs)
    en = end(exs)
if (viewtype == "exons") {
    ys = lapply(gr, function(x) pp(length(x)))
    yvs = unlist(ys) #1+(0:(nel-1))/nel
}
else if (viewtype == "transcripts") {
    tnms = exs$tx_id
    ft = factor(tnms)
    yvs = (as.numeric(ft)-1)/length(levels(ft))
}
else stop("viewtype not %in% c('exons', 'transcripts')")
    newdf = data.frame(st, en, yv = yvs, sym = sym)
    rng = range(exs)
    df = data.frame(range = c(start(rng), end(rng)), yval = rep(yval,2)) 
    strn = as.character(strand(exs)[1])
    ardir = ifelse(strn=="+", "last", "first")
    pl = ggplot(df, aes(x = range, y = yval)) + 
      geom_segment(aes(x = st, y = yv, xend = en, yend = yv, colour = sym),
           data = newdf, arrow=arrow(ends=ardir, length=unit(arrmm, "mm")), size=2.5)
    pl + xlab(as.character(GenomeInfoDb::seqnames(exs)[1]))
}

