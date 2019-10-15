from os import path
from os import environ
from glob import glob
import random
import string
import json
import pickle

from base64 import b64encode
import matplotlib.pyplot as plt


def random_string(n=5):
  return ''.join([random.choice(list(string.ascii_uppercase + string.digits)) for _ in range(n)])


class Granatum:

  def __init__(self, swd='/data'):
    self.swd = environ.get('GRANATUM_SWD', swd)
    self.uploaded_files_dir = path.join(self.swd, 'uploaded_files')
    self.exports_dir = path.join(self.swd, 'exports')
    self.exports_anno_file = path.join(self.swd, 'exports_anno.json')
    self.imports_dir = path.join(self.swd, 'imports')
    self.args_file = path.join(self.swd, 'args.json')
    self.results_file = path.join(self.swd, 'results.json')
    self.dynamic_exports = []
    self.results = []

    self.debug_dir = path.join(self.swd, 'debug')

    with open(path.join(self.swd, 'args.json')) as f:
      try:
        self.args = json.load(f)
      except FileNotFoundError:
        self.args = {}

  #-- uploaded_files -------------------------------------------------

  def get_uploaded_file_path(self, inject_into):
    try:
      return glob(path.join(self.uploaded_files_dir, inject_into, '*'))[0]
    except IndexError:
      return None

  #-- imports -------------------------------------------------

  def get_import(self, inject_into):
    import_file = path.join(self.imports_dir, inject_into)
    with open(import_file, 'r') as f:
      return json.load(f)

  #-- args -------------------------------------------------

  def get_arg(self, inject_into, default=None):
    return self.args.get(inject_into, default)

  #-- exports  -------------------------------------------------

  def export_statically(self, data, extract_from):
    with open(path.join(self.exports_dir, extract_from), 'w') as f:
      json.dump(data, f)

  def export(self, data, extract_from, kind=None, meta=None, dynamic=True):
    if dynamic:
      self.dynamic_exports.append({
        'extractFrom': extract_from,
        'kind': kind,
        'meta': meta,
      })

    self.export_statically(data, extract_from)


  #-- results  -------------------------------------------------

  def add_current_figure_to_results(self, description=None, zoom=2, width=750, height=650, dpi=100):
    save_filepath = path.join('/tmp', random_string() + '.png')

    fig = plt.gcf()
    fig.set_figheight(height / dpi)
    fig.set_figwidth(width / dpi)
    fig.savefig(save_filepath, dpi=zoom * dpi)

    with open(save_filepath, 'rb') as f:
      image_b64 = b64encode(f.read()).decode('utf-8')

    self.results.append(
      {
        'type': 'png',
        'width': width,
        'height': height,
        'description': description,
        'data': image_b64,
      }
    )

  def add_result(self, data, data_type='json', **kwargs):
    self.results.append({
      'type': data_type,
      'data': data,
      **kwargs,
    })

  #-- commit  -------------------------------------------------

  def commit(self):
    with open(self.exports_anno_file, 'w') as f:
      json.dump(self.dynamic_exports, f)

    with open(self.results_file, 'w') as f:
      json.dump(self.results, f)

  #-- error  -------------------------------------------------

  def error(self, message):
    raise Exception(message)

  #-- debug utils  -------------------------------------------------

  def _pickle(self, data, filename='debug.pickle'):
    with open(path.join(self.debug_dir, filename), 'wb') as f:
      pickle.dump(data, f)
