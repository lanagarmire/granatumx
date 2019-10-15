read_count <- function (raw_count = raw_count, out_dir = out_dir) 
{
    raw_count = as.matrix(raw_count)
    totalCounts_by_cell = colSums(raw_count)
    saveRDS(totalCounts_by_cell, file = paste0(out_dir, "totalCounts_by_cell.rds"))
    totalCounts_by_cell[totalCounts_by_cell == 0] = 1
    raw_count = sweep(raw_count, MARGIN = 2, totalCounts_by_cell/10^6, 
        FUN = "/")
    if (min(raw_count) < 0) {
        #print("smallest count cannot be negative!")
        stop()
    }
    count_lnorm = log10(raw_count + 1.01)
    return(count_lnorm)
}
