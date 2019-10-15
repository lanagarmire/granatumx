import Ajv from 'ajv';
import chalk from 'chalk';
import { readFileSync } from 'fs';
import yaml from 'js-yaml';
import Knex from 'knex';
import _ from 'lodash';

const ajv = new Ajv();
const gboxValidator = ajv.compile(yaml.safeLoad(readFileSync('./jsonSchemas/IGboxSpec.json', 'utf8')) as object);

export default async (knex: Knex, defaultMeta: object, gboxSpec: IGboxSpec) => {
  if (!gboxValidator(gboxSpec)) {
    console.log(chalk.bgRedBright.black('Error: Incorrect specification format!'));
    console.error(gboxValidator.errors);
    process.exit();
  }

  console.log(chalk.greenBright(`--->   Sending gbox ${gboxSpec.id} to database ...`));

  const combinedMeta = _.defaults(gboxSpec.meta, defaultMeta);

  const doesExist =
    +(await knex('gbox')
      .count()
      .where({ id: gboxSpec.id }))[0].count > 0;

  if (doesExist) {
    console.log(chalk.greenBright(`--->   Gbox ${gboxSpec.id} does exist in the database, updating ...`));

    await knex('gbox')
      .where({
        id: gboxSpec.id,
      })
      .update({
        meta: JSON.stringify(combinedMeta),
        endpoints: JSON.stringify(gboxSpec.endpoints),
        frontend: JSON.stringify(gboxSpec.frontend),
      });
  } else {
    console.log(chalk.greenBright(`--->   Gbox ${gboxSpec.id} does not exist in the database, writing ...`));

    await knex('gbox').insert({
      id: gboxSpec.id,
      meta: JSON.stringify(combinedMeta),
      endpoints: JSON.stringify(gboxSpec.endpoints),
      frontend: JSON.stringify(gboxSpec.frontend),
    });
  }
};
