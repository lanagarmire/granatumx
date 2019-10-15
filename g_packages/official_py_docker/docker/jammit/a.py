f = open('./__test/jammit.pickle', 'rb')

import pickle

jammit = pickle.load(f)

u = jammit.result_rows[11]['u']
print(u[u.nonzero()])
