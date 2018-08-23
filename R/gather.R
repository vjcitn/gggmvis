#' harvest encode690 components of AnnotationHub by ChIPped factor
#' @importFrom GenomeInfoDb seqlevelsStyle "seqlevelsStyle<-"
#' @param factor character(1) 
#' @param filtrange an optional GRanges that will be used to confine
#' information returned, using subsetByOverlaps
#' @note Uses encode690 to find AnnotationHub elements.  This returns
#' a list of GRanges as opposed to a GRangesList.  As of August 2018
#' there is a bug in GRangesList metadata() handling, so we work
#' with pure lists for now.
#' @return a list of GRanges, elements as provided by AnnotationHub
#' @examples
#' l1 = enc690ByFactor()
#' sapply(l1,length)
#' gm = GRanges("chr17", IRanges(38077296, 38083884))
#' \dontrun{
#' l2 = lapply(l1, function(x) subsetByOverlaps(x, gm+10000))
#' l3 = enc690ByFactor(filtrange=gm+10000)
#' }  # that checks out
#' # here we demonstrate how to use output to obtain
#' # material for visualization of cell-type specific
#' # binding near a given gene
#' br = GRanges("chr13", IRanges(32889617, 32973809)) # BRCA2
#' # could do following but slow ...
#' # dd = enc690ByFactor(filtrange=br+10000)
#' dd = lapply(l1, function(x) subsetByOverlaps(x, br+10000))
#' lens = sapply(dd,length)
#' cls = sapply(dd, function(x) metadata(x)$cell)
#' cls = rep(cls, lens)
#' ee = do.call(rbind, lapply(dd, as.data.frame))
#' ee$cell = factor(as.character(cls))
#' ee$yval = 1+(as.numeric(factor(as.character(cls)))-1)/length(unique(cls))
#' edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
#' ggvisForSymbol("BRCA2", resource=edb) +
#'    geom_segment(aes(x=start, xend=end, y=yval, yend=yval,
#'          group=cell, colour=cell), data=ee, size=2.5) +
#'    theme(axis.text.y = element_blank(), axis.title.y=element_blank()) + 
#'        ylim(-.5,2) + ggtitle("CEBPB binding near BRCA2")
#' @export
enc690ByFactor = function (factor = "CEBPB", filtrange=NULL) 
{
    data(encode690)
    encode690 = as.data.frame(encode690)
    stopifnot(factor %in% encode690$target)
    tmp = dplyr::filter(encode690, target == factor)[, c("AHID", "cell")]
    ids = tmp[["AHID"]]
    cells = tmp[["cell"]]
  message(paste("retrieving", length(ids), "AnnotationHub elements."))
  suppressWarnings({
   suppressMessages({
    ans = lapply(seq_len(length(ids)), function(x) {
        cat(".")
        tmp = AnnotationHub::AnnotationHub()[[ids[x]]]
        if (!is.null(filtrange)) {
          seqlevelsStyle(filtrange) = seqlevelsStyle(tmp)
          tmp = subsetByOverlaps(tmp, filtrange)
          }
        metadata(tmp) = c(metadata(tmp), cell = cells[x])
        tmp
    })
    })
    })
    cat("\n")
    names(ans) = ids
    ans
}
