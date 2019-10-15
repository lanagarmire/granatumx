write_count <-
function (count_imp, out_dir = out_dir) 
{
    totalCounts_by_cell = readRDS(paste0(out_dir, "totalCounts_by_cell.rds"))
    count_imp = sweep(count_imp, MARGIN = 2, totalCounts_by_cell/10^6, 
        FUN = "*")
    count_imp = round(count_imp, digits = 2)

    return(count_imp)
}
