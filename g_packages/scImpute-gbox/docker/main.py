#!/usr/bin/env python3

from sys import stdin,stdout,exit
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import json
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from base64 import b64encode
from glob import glob
from granatum_sdk import Granatum

def main():
    gn = Granatum()

gn.get_import("assay")
    input_args = {
        'raw_count': np.array(input_assay.get('matrix')),
        'drop_thre': gn.get_arg('drop_thre')
        }

    rpy2.rinterface.set_writeconsole_regular(lambda x: None)
    rpy2.rinterface.set_writeconsole_warnerror(lambda x: None)
    rpy2.robjects.numpy2ri.activate()
    for r_file in glob('scImpute/*.R'):
        robjects.r('source("{}")'.format(r_file))

    r_scimpute = robjects.globalenv['scimpute']
    imputated_array = np.array(r_scimpute(**input_args))

    # Plot raw and imputated matrices
    fig, ax = plt.subplots(1,2)
    sns.heatmap(np.log(1+input_args['raw_count']), ax=ax[0], xticklabels=False, yticklabels=False)
    sns.heatmap(np.log(1+imputated_array), ax=ax[1], xticklabels=False, yticklabels=False)
    ax[0].set_title('Raw count',fontweight='bold',fontsize=14)
    ax[0].set_xlabel('Genes',fontsize=12)
    ax[0].set_ylabel('Cells',fontsize=12)
    ax[1].set_title('Imputated counts (log10)',fontweight='bold',fontsize=14)
    ax[1].set_xlabel('Genes',fontsize=12)
    ax[1].set_ylabel('Cells',fontsize=12)

    gn.add_current_figure_to_results(
        "Expression levels comparison: before vs. after imputation",
        zoom=3,
        dpi=50,
        height=400,
    )

    gn.add_result(
        "\n".join(
            [
                "Percentage of zeros before filtering: **{}**".format(np.sum(imputated_array==0,axis=0)/imputated_array.shape[0]),
                "",
                "Percentage of zeros after filtering: **{}**".format(np.sum(input_args['raw_count']==0,axis=0)/input_args['raw_count'].shape[0]),
            ]
        ),
        type="markdown",
    )

    input_assay['matrix'] = imputated_array.tolist()
    gn.export(input_assay, "Inputed Assay", dynamic=False)

    gn.commit()


if __name__ == '__main__':
    main()
