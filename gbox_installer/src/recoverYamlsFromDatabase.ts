import fs from 'fs';
import { ensureDirSync } from 'fs-extra';
import glob from 'glob';
import yaml from 'js-yaml';
import Knex from 'knex';

ensureDirSync('recovered');

(async () => {
  const knex = await Knex({
    client: 'pg',
    connection: process.env.DATABASE_URL || 'postgres://postgres:12qw@localhost:5433/granatum',
  });

  const gboxes = await knex('gbox').select(['id', 'meta', 'endpoints', 'frontend']);

  await Promise.all(
    gboxes.map(async (gbox: any) => {
      fs.writeFileSync(
        `recovered/${gbox.id}.yaml`,
        yaml.safeDump(
          JSON.parse(
            JSON.stringify({
              meta: gbox.meta,
              endpoints: gbox.endpoints,
              frontend: gbox.frontend,
            }),
          ),
        ),
      );
    }),
  );

  knex.destroy();
})();
