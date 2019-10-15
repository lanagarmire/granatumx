import chalk from 'chalk';
import { readFileSync } from 'fs';
import glob from 'glob';
import yaml from 'js-yaml';
import Knex from 'knex';
import * as path from 'path';
import installOnePackage from './installOnePackage';
import installOneRecipe from './installOneRecipe';

const listOfPackagesPath = './listOfPackages.yaml';
const recipeYamls = glob.sync('./recipes/*.yaml');

(async () => {
  const knex = await Knex({
    client: 'pg',
    connection: process.env.DATABASE_URL || 'postgres://postgres:12qw@localhost:5433/granatum',
  });

  console.log(chalk.blueBright(`===> Install packages`));
  const listOfPackages = yaml.safeLoad(readFileSync(listOfPackagesPath, 'utf8')) as string[];
  for (const p of listOfPackages) {
    const fullPath = path.resolve(p);
    const packageSpecPath = path.resolve(fullPath, 'package.yaml');
    console.log(chalk.blueBright(`===>   Loading package spec from ${packageSpecPath}`));
    const packageSpec = yaml.safeLoad(readFileSync(packageSpecPath, 'utf8')) as IPackageSpec;
    await installOnePackage(knex, fullPath, packageSpec, { verbose: true });
  }

  console.log(chalk.blueBright(`===> Install recipes`));
  await knex('recipe_gbox').delete();
  for (const y of recipeYamls) {
    console.log(chalk.blueBright(`===>   Loading recipe spec from ${y}`));
    const recipeSpec = yaml.safeLoad(readFileSync(y, 'utf8')) as IRecipeSpec;
    await installOneRecipe(knex, recipeSpec);
  }

  knex.destroy();

  console.log(chalk.blueBright(`.... All done!`));
})();
