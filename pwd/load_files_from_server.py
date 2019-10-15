#!/usr/bin/env python

from os.path import basename

import pandas as pd

from granatum_sdk import Granatum

#import numpy as np

def main():
  gn = Granatum()
  assay_file = gn.get_arg('pathToAssayFile')
  sample_meta_file = gn.get_arg('pathToSampleMetaFile')

  tb = pd.read_csv(assay_file, index_col=0)
  sample_ids = tb.columns.values.tolist()

  assay_export_name = '[A]{}'.format(basename(assay_file))

  exported_assay = {
    'matrix': tb.values.tolist(),
    'sampleIds': sample_ids,
    'geneIds': tb.index.values.tolist(),
  }

  gn.export(exported_assay, assay_export_name, 'assay')

  results_markdown = 'The assay has **{}** genes and **{}** samples.'.format(tb.shape[0], tb.shape[1])
  results_markdown += '\n\n The first few rows and columns:'
  results_markdown += '\n\n```'
  results_markdown += '\n'.join([', '.join(x) for x in tb.values[:10, :10].astype(str).tolist()])
  results_markdown += '\n```'

  if sample_meta_file is not None:
    results_markdown += '\n'
    sample_meta_tb = pd.read_csv(sample_meta_file)

    for meta_name in sample_meta_tb.columns:
      meta_output_name = '[M]{}'.format(meta_name)

      sample_meta_dict = dict(zip(sample_ids, sample_meta_tb[meta_name].values.tolist()))

      gn.export(sample_meta_dict, meta_output_name, 'sampleMeta')

      num_sample_values = 5
      sample_values = ', '.join(sample_meta_tb[meta_name].astype(str).values[0:num_sample_values].tolist())
      num_omitted_values = len(sample_meta_tb[meta_name]) - num_sample_values

      if num_omitted_values > 0:
        etc = ', ... and {} more entries'.format(num_omitted_values)
      else:
        etc = ''

      results_markdown += '\n  * Sample meta with name **{}** is accepted ({}{}).'.format(meta_name, sample_values, etc)

  gn.add_result(results_markdown, 'markdown')

  #TODO: SAVE assay pickle

  gn.commit()


if __name__ == '__main__':
  main()
