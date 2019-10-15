#!/usr/bin/env python3

print("""
Hello, welcome to Cedric's scImpute docker image.

Available commands:

  scImpute - docker run -i cedric-scimpute python3 ./main.py

               stdin format: {
                 assay: [matrix] (required),
                 drop_thre: [number between 0 and 1] (default:0.5)
               }
""")
