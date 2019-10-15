import numpy as np

from granatum_sdk import Granatum


def main():
    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import("assay"))

    frob_norm = np.linalg.norm(df.values)

    df = df / frob_norm

    gn.add_result(
        f"""\
The original assay had Frobenius norm of {frob_norm}, after normalization its
Frobenius norm is now {np.linalg.norm(df.values)}""",
        'markdown',
    )

    gn.export(gn.assay_from_pandas(df), "Frobenius normalized assay", dynamic=False)

    gn.commit()


if __name__ == "__main__":
    main()
