import time

import numpy as np
import pandas as pd
import random

from .ssvd import ssvd
from .utils import array_visualizer, center_rows, eprint, permute_each_row
from multiprocessing import Pool


def split_array_by_track(v, track_vec):
    return dict([(k, v.values) for k, v in pd.Series(v).groupby(track_vec)])


def do_for_one_permute(task):
    np.random.seed()
    thetas = task['thetas']
    xx = task['xx']
    feature_meta = task['feature_meta']
    verbose = task['verbose']
    kwargs = task['kwargs']
    df_max = np.abs(xx).max()

    permuted_xx = permute_each_row(xx)
    tu, _, _ = np.linalg.svd(permuted_xx, full_matrices=False)

    n_sigs_list = []
    n_sigs_by_track_list = []
    for i_theta, theta in enumerate(thetas):
        alpha = df_max * theta,
        # nfps (numbers of false positives) is a list of the numbers of
        # false positive features in the permuted runs Any positive given
        # by ssvd in a permuted matrix is deem a false positive
        fu, _ = ssvd(permuted_xx, tu=tu, alpha=alpha, verbose=verbose, **kwargs)
        n_sigs_list.append(np.count_nonzero(fu))
        fu_by_track = split_array_by_track(fu, feature_meta['track'].values)
        n_sigs_by_track_list.append({track_name: np.count_nonzero(fu) for track_name, fu in fu_by_track.items()})

    n_sigs_vec = np.array(n_sigs_list)
    n_sigs_vec_by_track = {k: np.array([d[k] for d in n_sigs_by_track_list]) for k in feature_meta['track'].unique()}

    return {
        'id_': 'permute',
        'n_sigs_vec': n_sigs_vec,
        'n_sigs_vec_by_track': n_sigs_vec_by_track,
    }


class JAMMIT:

    def __init__(self, df, feature_meta, verbose=1):
        self.verbose = verbose

        # store the number of rows (features) in each matrix
        self.feature_meta = feature_meta

        # stack matrices
        self.df = df

        self.result_rows = None
        self.performance_stats = {}

    def center_data(self):
        self.df.values = center_rows(self.df.values)

    @classmethod
    def from_df_and_track_vec(cls, df, track_vec, **kwargs):
        return cls(df, pd.DataFrame({'track': track_vec}))

    # `df_dict` initializer (df_dict is a dict df_name -> df)
    @classmethod
    def from_df_dict(cls, df_dict, **kwargs):
        # store the number of rows (features) in each matrix
        feature_meta = pd.concat([pd.DataFrame({'track': df_name}, index=df.index) for df_name, df in df_dict.items()])

        # stack matrices
        df = pd.concat(list(df_dict.values()))

        return cls(df, feature_meta)

    @classmethod
    def from_csvs(
        cls,
        paths,
        sep=",",
        header=None,
        index_col=None,
        genes_as_columns=False,
        **kwargs,
    ):
        df_dict = {}
        for track_name, p in paths.items():
            df = pd.read_csv(p, sep=sep, header=header, index_col=index_col)

            if genes_as_columns:
                df = df.transpose()

            df_dict[track_name] = df

        return cls.from_df_dict(df_dict, **kwargs)

    @classmethod
    def from_dfs(cls, dfs, **kwargs):
        return cls.from_df_dict({f"matrix_{i_df}": df for i_df, df in enumerate(dfs)}, **kwargs)

    def get_alpha_from_theta(self, theta):
        return np.abs(self.df.values).max() * theta

    def run_for_one_alpha(self, alpha, verbose=None, **kwargs):
        tik = time.time()

        verbose = self.verbose if verbose is None else verbose

        u, v = ssvd(self.df.values, alpha=alpha, verbose=verbose, **kwargs)

        return {
            'u': u,
            'v': v,
        }

    def scan(
        self,
        thetas,
        calculate_fdr=True,
        n_perms=2,
        verbose=None,
        **kwargs,
    ):
        tik = time.time()

        verbose = self.verbose if verbose is None else verbose

        result_rows = []
        for i_theta, theta in enumerate(thetas):
            if verbose >= 1:
                eprint(f'({i_theta}/{len(thetas)}) theta = {theta}')

            alpha = self.get_alpha_from_theta(theta)

            u, v = ssvd(self.df.values, alpha=alpha, verbose=verbose, **kwargs)
            n_sigs = np.count_nonzero(u)
            n_sigs_by_track = {
                track_name: np.count_nonzero(uu)
                for track_name, uu in split_array_by_track(u, self.feature_meta['track'].values).items()
            }


            result_rows.append(
                {
                    "theta": theta,
                    "alpha": alpha,
                    "u": u,
                    "v": v,
                    "n_sigs": n_sigs,
                    'n_sigs_by_track': n_sigs_by_track,
                }
            )

        if calculate_fdr:
            with Pool(1) as p:
                do_for_one_permute_list = []
                for i_perm in range(n_perms):
                    do_for_one_permute_list.append(
                        {
                            'xx': self.df.values,
                            'feature_meta': self.feature_meta,
                            'thetas': thetas,
                            'verbose': verbose,
                            'kwargs': kwargs,
                        }
                    )

                results_iter = p.imap(do_for_one_permute, do_for_one_permute_list)
                n_sigs_vec_list = []
                n_sigs_vec_by_track_list = []
                for i_result, result in enumerate(results_iter):
                    if verbose >= 1:
                        eprint(f"({i_result}/{n_perms}) {result['id_']}")
                    n_sigs_vec_list.append(result['n_sigs_vec'])
                    n_sigs_vec_by_track_list.append(result['n_sigs_vec_by_track'])

            n_sigs_mat = np.array(n_sigs_vec_list)
            n_sigs_mat_by_track = {
                track_name: np.array([d[track_name] for d in n_sigs_vec_by_track_list
                                     ]) for track_name in self.feature_meta['track'].unique()
            }

            expected_nfp_by_theta = np.mean(n_sigs_mat, axis=0)
            expected_nfp_by_theta_by_track = {
                track_name: np.mean(n_sigs_mat, axis=0) for track_name, n_sigs_mat in n_sigs_mat_by_track.items()
            }

            for i_theta, _ in enumerate(thetas):
                fdr = np.minimum(expected_nfp_by_theta[i_theta] / n_sigs, 1)
                fdr_by_track = {
                    track_name:
                    np.minimum(expected_nfp_by_theta_[i_theta] / n_sigs_by_track[track_name], 1)
                    for track_name, expected_nfp_by_theta_ in expected_nfp_by_theta_by_track.items()
                }
                result_rows[i_theta]['fdr'] = fdr
                result_rows[i_theta]['fdr_by_track'] = fdr_by_track

        tok = time.time()
        time_elapsed = tok - tik

        self.performance_stats = {
            "execution_duration": time_elapsed,
        }
        self.result_rows = result_rows

    def format(
        self, columns=None, show_vector_preview=False, show_individual_tracks=False, to_csv_path="./result_table.csv"
    ):
        result_table = pd.DataFrame(self.result_rows)

        formatted_table = pd.DataFrame()

        formatted_table["theta"] = result_table["theta"]
        formatted_table["alpha"] = result_table["alpha"]
        formatted_table["n_sigs"] = result_table["n_sigs"]
        formatted_table["fdr"] = result_table["fdr"]

        if show_vector_preview:
            formatted_table["u"] = result_table["u"].map(array_visualizer)
            formatted_table["v"] = result_table["v"].map(array_visualizer)

        if show_individual_tracks:
            formatted_table['n_sigs_by_track'] = [str(x) for x in result_table["n_sigs_by_track"]]
            formatted_table['fdr_by_track'] = [str(x) for x in result_table["fdr_by_track"]]
            # for track_name in self.feature_meta['track']:
            #     formatted_table[f"n_sigs_of_{track_name}"] = [d[track_name] for d in result_table["n_sigs_by_track"]]
            #     formatted_table[f"fdr_of_{track_name}"] = [d[track_name] for d in result_table["fdr_by_track"]]

        if type(to_csv_path) is str:
            formatted_table.to_csv(to_csv_path)

        if columns is None:
            return formatted_table
        else:
            return formatted_table[columns]
