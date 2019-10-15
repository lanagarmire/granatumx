main:
	docker-compose up -d && sleep 1 && make db-setup

db-setup:
	cat setup.sql | docker run -i --rm --network=host postgres:10.3 psql -h localhost -p 5433 -U postgres -d postgres -f -
