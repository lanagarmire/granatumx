import fs from 'fs';
import glob from 'glob';
import yaml from 'js-yaml';
import Knex from 'knex';

const gboxYamls = process.argv.slice(2);

(async () => {
  const knex = await Knex({
    client: 'pg',
    connection: process.env.DATABASE_URL || 'postgres://postgres:12qw@localhost:5433/granatum',
  });

  console.log(`==> updaing ${gboxYamls.length} gbox(es):`);

  await Promise.all(
    gboxYamls.map(async (y) => {
      try {
        const gboxName = y.replace(/^.*\/([^/]+)\.yaml$/, '$1');
        const spec:
          | { meta?: object; endpoints?: Array<string | object>; frontend?: string | object }
          | undefined = yaml.safeLoad(fs.readFileSync(y, 'utf8'));

        if (spec == null) {
          // noinspection ExceptionCaughtLocallyJS
          throw new Error('The spec is null.');
        }

        const doesExist =
          +(await knex('gbox')
            .count()
            .where({ id: gboxName }))[0].count > 0;

        if (doesExist) {
          console.log(`${gboxName} (overwrite)`);

          await knex('gbox')
            .where({
              id: gboxName,
            })
            .update({
              meta: JSON.stringify(spec.meta),
              endpoints: JSON.stringify(spec.endpoints),
              frontend: JSON.stringify(spec.frontend),
            });
        } else {
          console.log(`${gboxName} (new)`);

          await knex('gbox').insert({
            id: gboxName,
            meta: JSON.stringify(spec.meta),
            endpoints: JSON.stringify(spec.endpoints),
            frontend: JSON.stringify(spec.frontend),
          });
        }
      } catch (e) {
        console.log(e);
      }
    }),
  );

  knex.destroy();
})();
