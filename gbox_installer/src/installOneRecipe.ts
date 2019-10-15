import Ajv from 'ajv';
import chalk from 'chalk';
import { readFileSync } from 'fs';
import yaml from 'js-yaml';
import Knex from 'knex';

const ajv = new Ajv();
const recipeValidator = ajv.compile(yaml.safeLoad(readFileSync('./jsonSchemas/IRecipeSpec.json', 'utf8')) as object);

export default async (knex: Knex, recipeSpec: IRecipeSpec) => {
  if (!recipeValidator(recipeSpec)) {
    console.log(chalk.bgRedBright.black('Error: Incorrect specification format!'));
    console.error(recipeValidator.errors);
    process.exit();
  }

  const doesExist =
    +(await knex('recipe')
      .count()
      .where({ id: recipeSpec.id }))[0].count > 0;

  if (doesExist) {
    console.log(chalk.greenBright(`---> ${recipeSpec.meta.title} (overwrite)`));

    await knex('recipe')
      .where({
        id: recipeSpec.id,
      })
      .update({
        meta: JSON.stringify(recipeSpec.meta),
      });
  } else {
    console.log(chalk.greenBright(`---> ${recipeSpec.meta.title} (new)`));

    await knex('recipe').insert({
      id: recipeSpec.id,
      meta: JSON.stringify(recipeSpec.meta),
    });
  }

  if (recipeSpec.steps != null) {
    await Promise.all(
      recipeSpec.steps.map(async (step, i) => {
        if (recipeSpec.meta != null) {
          console.log(chalk.greenBright(`--->   - ${step.gbox}`));
        }

        await knex('recipe_gbox').insert({
          recipe_id: recipeSpec.id,
          gbox_id: step.gbox,
          rank: i,
          initial_state: step.initialState,
        });
      }),
    );
  }
};
