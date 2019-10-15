main: build

build:
	docker build -t granatum/lncrna .

shell:
	docker run -it granatum/lncrna /bin/bash

exec: build
	docker run -it granatum/lncrna Rscript ./test.R
