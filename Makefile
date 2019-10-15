# prerequisites: docker yarn node

setup: stop webapp-setup gboxes-setup task-runner-setup folder-setup db-init start-dev

update: stop webapp-setup gboxes-setup task-runner-setup folder-setup start-dev

db-up:
	docker container inspect temporary-postgres-instance > /dev/null 2>&1 \
    || docker run -it --rm -v /var/granatum/postgres_data:/var/lib/postgresql/data:rw -p 5433:5432 -e POSTGRES_PASSWORD=12qw --name temporary-postgres-instance -d postgres:10.9 -c 'listen_addresses=*'

start: db-up
	tmux new-session -s GranatumX 'make webapp-prod; bash' \; split-window -h 'make task-runner; bash'

start-dev: db-up
	tmux new-session -s GranatumX 'make webapp; bash' \; split-window -h 'make task-runner; bash'

stop: db-up
	tmux kill-session -t GranatumX; true

# for production: remove make db-mock
db-init: db-up
	sleep 4 && make db-setup \
					&& sleep 1 && cd gbox_installer && node ./dist/index.js \
					&& sleep 1 && cd .. && make db-mock

folder-setup:
	mkdir -p /var/granatum

task-runner:
	cd ./taskRunner && yarn start

webapp:
	cd ./webApp && yarn start

webapp-prod:
	cd ./webApp/build-prod && env NODE_ENV=production yarn start


webapp-setup:
	cd ./webApp && yarn install

gboxes-setup:
	cd ./gbox_installer && yarn install

task-runner-setup:
	cd ./taskRunner && yarn install

psql:
	docker run -it --rm --network=host postgres:10.9 psql -h localhost -p 5433 -U postgres -d granatum

db-setup:
	cat db/setup.sql | docker run -i --rm -e PGPASSWORD=12qw --network=host postgres:10.9 psql -h localhost -p 5433 -U postgres -d postgres -f -

db-mock:
	cat db/data.sql | docker run -i --rm -e PGPASSWORD=12qw --network=host postgres:10.9 psql -h localhost -p 5433 -U postgres -d granatum -f -

# old:
# 	cat setup.sql data.sql | docker run -i --rm --link granatum-postgres-instance:postgres -e PGPASSWORD=12qw granatum/postgres psql -h postgres -U postgres -d postgres -f - \
#             && postgraphql -c postgres://postgres:12qw@localhost:5433/granatum_mock_2 -s public -p 3000 -r anonymous -e ufdMnOWZ41 \
#                            -t public.jwt_token --show-error-stack --extended-errors hint,detail,errcode --watch

postgraphile:
	postgraphile -c postgres://postgres:12qw@localhost:5433/granatum -s public -p 3000 -r anonymous -e ufdMnOWZ41 \
                           -t public.jwt_token --disable-default-mutations --show-error-stack --extended-errors hint,detail,errcode --watch

update-submodules:
	git submodule update --init --recursive

dockerd-start:
	@echo -e "\033[0;33m##### Make sure to run with sudo #####\033[0m"
	tmux new-session -s dockerd 'dockerd'

status:
	cd ./webApp && yarn status

itest:
	cd ./webApp && yarn itest
