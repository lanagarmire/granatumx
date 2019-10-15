import json
from sys import stdin, stdout

# Retrieve and parse the JSON from stdin
input = json.load(stdin)

# Disassemble the input for function invocation
args = {
  'arg1': input.get('arg1'),
  'arg2': input.get('arg2'),
  ...
}

# Try running the main function of your package
try:
  results = your_function(**args)
  # If the run is successful, return the results
  output = {
    exports: {...},
    results: {...}
  }
except Error:
  # If there were errors, report them to Granatum
  output = {
    errors: [...],
  }

# Compose and write the output JSON into stdout
json.dump(output, stdout)
