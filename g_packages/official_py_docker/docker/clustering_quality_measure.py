#!/usr/bin/env python

import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.metrics import adjusted_rand_score

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    sample_meta_true = gn.get_import("sample_meta_true")
    sample_meta_predicted = gn.get_import("sample_meta_predicted")

    # Using pandas series to align the two metas in case they have different sample IDs
    rand_score = adjusted_rand_score(
        pd.Series(sample_meta_true), pd.Series(sample_meta_predicted)
    )
    mutual_info_score = adjusted_mutual_info_score(
        pd.Series(sample_meta_true), pd.Series(sample_meta_predicted)
    )

    results_markdown = "\n".join(
        [
            "Adjusted Rand score: **{}**".format(rand_score),
            "",
            "Adjusted mutual information score: **{}**".format(mutual_info_score),
        ]
    )

    gn.add_result(results_markdown, "markdown")
    gn.commit()


if __name__ == "__main__":
    main()
