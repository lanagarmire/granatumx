library(tidyverse)
library(gbm)

source('./granatum_sdk.R')

GB_model_1 <- readRDS('./GB_model_1.RDS')

n_trees <- gn_get_arg('n_trees')

pred_df <- data_frame(
  TMPO.AS1.202    = gn_get_arg('TMPO.AS1.202'),
  SBF2.AS1.201    = gn_get_arg('SBF2.AS1.201'),
  MAFG.AS1.201    = gn_get_arg('MAFG.AS1.201'),
  LINC01852.201   = gn_get_arg('LINC01852.201'),
  LINC01614.201   = gn_get_arg('LINC01614.201'),
  LINC00261.202   = gn_get_arg('LINC00261.202'),
  HHIP.AS1.203    = gn_get_arg('HHIP.AS1.203'),
  HHIP.AS1.201    = gn_get_arg('HHIP.AS1.201'),
  CARD8.AS1.201   = gn_get_arg('CARD8.AS1.201'),
  AL109741.1.201  = gn_get_arg('AL109741.1.201'),
  AC146944.4.201  = gn_get_arg('AC146944.4.201'),
  AC027288.3.203  = gn_get_arg('AC027288.3.203'),
  LINC02555.201   = gn_get_arg('LINC02555.201'),
  LINC01936.201   = gn_get_arg('LINC01936.201'),
  TBX5.AS1.201    = gn_get_arg('TBX5.AS1.201'),
  HSPC324.201     = gn_get_arg('HSPC324.201'),
  GATA6.AS1.202   = gn_get_arg('GATA6.AS1.202'),
  AP000866.2.201  = gn_get_arg('AP000866.2.201'),
  ADAMTS9.AS2.203 = gn_get_arg('ADAMTS9.AS2.203'),
  AC008268.1.201  = gn_get_arg('AC008268.1.201')
)

GB_model_1 <- readRDS('./GB_model_1.RDS')
GB_val_1 <- predict.gbm(GB_model_1, pred_df, n.trees=n_trees, type='response', single.tree = TRUE)
## gn_add_result(sprintf('GB_val_1 = %s', GB_val_1), 'markdown', list(label='Debug', description='Method selected to normalize the assay.'))

p <- ggplot() +
    geom_col(aes(x='', y=c(1 - GB_val_1, GB_val_1), fill=fct_inorder(c('No risk', 'Lung cancer risk'))), width=1) +
    geom_text(aes(x='', y=c(GB_val_1 / 2, 1 - (1 - GB_val_1)/2), label=c(sprintf('Lung cancer risk (%s)', scales::percent_format()(GB_val_1)), sprintf('No risk (%s)', scales::percent_format()(1 - GB_val_1)))), width=1, color='white', size=10) +
    scale_y_continuous() +
    scale_fill_manual(values=c('lightgrey', 'red'), guide=F) +
    coord_polar(theta="y", start=0, direction=1) +
    theme_void()
gn_add_ggplot_to_results(p, 'Percentage of lung cancer risk', height=700, dpi=100)

gn_commit()
