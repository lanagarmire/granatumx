from granatum_clustering.clustering_apps import GranatumClustering

from sys import stdin
from sys import stdout

import json

import numpy as np


def main():
    inputs = json.load(stdin)

    granatum_clustering = GranatumClustering(**inputs)

    matrix = inputs['assay']['matrix']
    matrix = np.array(matrix).T

    sample_ids = inputs['assay']['sampleIds']

    results = granatum_clustering.fit(
        matrix=matrix, sample_ids=sample_ids, jsonify=True)

    stdout.write(json.dumps(results))


if __name__ == '__main__':
    main()
