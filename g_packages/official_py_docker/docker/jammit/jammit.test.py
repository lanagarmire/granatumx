import pickle

import numpy as np
from jammit import JAMMIT

jammit = JAMMIT.from_csvs({'main': './__test/datamatcfro.csv', 'superviser': './__test/stepcfro.csv'})

jammit.scan(
    thetas=np.linspace(0, 1, 21),
    calculate_fdr=True,
    n_perms=10,
    verbose=1,
    convergence_threshold=0.000000001,
)

print(jammit.format(to_csv_path="./__test/result_table.csv"))

for k, v in jammit.performance_stats.items():
    print(f'{k} = {v}')

with open("./__test/jammit.pickle", "wb") as f:
    pickle.dump(jammit, f)
