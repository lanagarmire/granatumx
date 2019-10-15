import Knex from 'knex';
import config from './config';

let knex;

export default async () => {
  if (knex === undefined) {
    knex = await Knex({
      client: 'pg',
      connection: config.databaseUrl,
    });
  }
  return knex;
};
