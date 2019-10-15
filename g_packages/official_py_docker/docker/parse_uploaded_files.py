#!/usr/bin/env python

from os.path import basename

import pandas as pd

from collections import Counter
from granatum_sdk import Granatum, guess_gene_id_type, biomart_col_dict, convert_gene_ids

import numpy as np


def main():
    gn = Granatum()

    assay_file = gn.get_uploaded_file_path("assayFile")
    sample_meta_file = gn.get_uploaded_file_path("sampleMetaFile")
    file_format = gn.get_arg("fileFormat")
    species = gn.get_arg("species")

    if file_format == "csv":
        tb = pd.read_csv(assay_file, sep=",", index_col=0)
    elif file_format == "tsv":
        tb = pd.read_csv(assay_file, sep="\t", index_col=0)
    elif file_format == "excel":
        tb = pd.read_excel(assay_file, index_col=0)
    else:
        gn.error("Unknown file format: {}".format(file_format))

    sample_ids = tb.columns.values.tolist()
    gene_ids = tb.index.values.tolist()

    gene_id_type = guess_gene_id_type(gene_ids[:5])

    whether_convert_id = gn.get_arg("whether_convert_id")

    if whether_convert_id:
        to_id_type = gn.get_arg("to_id_type")
        add_info = gn.get_arg("add_info")

        # if there are duplicated ids, pick the first row
        # TODO: Need to have a more sophisticated handling of duplicated ids

        gene_ids, new_meta = convert_gene_ids(gene_ids, gene_id_type, to_id_type, species, return_new_meta=True)

        # TODO: remove NaN rows
        # TODO: combine duplicated rows

        if add_info:
            for col_name, col in new_meta.iteritems():
                gn.export(col.to_dict(), col_name, "geneMeta")

    assay_export_name = "[A]{}".format(basename(assay_file))

    exported_assay = {
        "matrix": tb.values.tolist(),
        "sampleIds": sample_ids,
        "geneIds": gene_ids,
    }

    gn.export(exported_assay, assay_export_name, "assay")

    entry_preview = '\n'.join([', '.join(x) for x in tb.values[:10, :10].astype(str).tolist()])

    gn.add_result(
        f"""\
The assay has **{tb.shape[0]}** genes (with inferred ID type: {biomart_col_dict[gene_id_type]}) and **{tb.shape[1]}** samples.

The first few rows and columns:

```
{entry_preview}
```
""",
        "markdown",
    )

    meta_rows = []
    if sample_meta_file is not None:
        if file_format == "csv":
            sample_meta_tb = pd.read_csv(sample_meta_file)
        elif file_format == "tsv":
            sample_meta_tb = pd.read_csv(sample_meta_file, sep="\t")
        elif file_format == "excel":
            sample_meta_tb = pd.read_excel(sample_meta_file)
        else:
            gn.error("Unknown file format: {}".format(file_format))

        for meta_name in sample_meta_tb.columns:
            meta_output_name = "[M]{}".format(meta_name)

            sample_meta_dict = dict(zip(sample_ids, sample_meta_tb[meta_name].values.tolist()))

            gn.export(sample_meta_dict, meta_output_name, "sampleMeta")

            num_sample_values = 5
            sample_values = ", ".join(sample_meta_tb[meta_name].astype(str).values[0:num_sample_values].tolist())
            num_omitted_values = len(sample_meta_tb[meta_name]) - num_sample_values

            if num_omitted_values > 0:
                etc = ", ... and {} more entries".format(num_omitted_values)
            else:
                etc = ""

            meta_rows.append({
                'meta_name': meta_name,
                'sample_values': str(sample_values) + etc,
            })

    # meta_message = '\n'.join(
    #     "* Sample meta with name **{meta_name}** is accepted ({sample_values}).".format(**x) for x in meta_rows
    # )

    # gn.add_result(meta_message, "markdown")

    # gn.add_result({'columns': []}, 'table')

    # TODO: SAVE assay pickle

    gn.commit()


if __name__ == "__main__":
    main()
