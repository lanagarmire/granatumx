import Pool from 'pg-pool';
import { URL } from 'url';
import config from './config';

const params = new URL(config.databaseUrl);

const pgConfig = {
  host: params.hostname,
  user: params.username,
  password: params.password,
  port: +params.port,
  database: params.pathname.split('/')[1],
};

export default new Pool(pgConfig);
