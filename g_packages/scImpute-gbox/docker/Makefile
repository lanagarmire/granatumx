main:
	docker build -t cedric-scimpute .

run:
	docker run cedric-scimpute

shell:
	docker run -it cedric-scimpute /bin/bash

test:
	cat input.json | docker run -i cedric-scimpute python3 ./main.py
