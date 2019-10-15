import Ajv from 'ajv';
import chalk from 'chalk';
import { execSync } from 'child_process';
import { readFileSync } from 'fs';
import yaml from 'js-yaml';
import Knex from 'knex';
import * as path from 'path';
import installOneGbox from './installOneGbox';

const ajv = new Ajv();
const packageValidator = ajv.compile(yaml.safeLoad(readFileSync('./jsonSchemas/IPackageSpec.json', 'utf8')) as object);

export default async (
  knex: Knex,
  cwd: string,
  packageSpec: IPackageSpec,
  options: { verbose: boolean } = { verbose: false },
) => {
  if (!packageValidator(packageSpec)) {
    console.log(chalk.bgRedBright.black('Error: Incorrect specification format!'));
    console.error(packageValidator.errors);
    process.exit();
  }

  const { buildCommand, gboxes, meta, id } = packageSpec;

  console.log(chalk.greenBright(`---> Building package: ${id}`));

  if (buildCommand) {
    console.log(chalk.greenBright(`--->   Running the buildCommand ...`));
    execSync(buildCommand, { stdio: options.verbose ? 'inherit' : 'ignore', cwd });
  } else {
    console.log(chalk.yellowBright(`--->   No buildCommand found`));
  }

  await Promise.all(
    gboxes.map(async (gbox) => {
      const gboxSpec = (() => {
        if (typeof gbox === 'string') {
          const gboxFullPath = path.resolve(cwd, gbox);
          console.log(chalk.greenBright(`--->   Loading gbox spec at ${gboxFullPath} ...`));
          return yaml.safeLoad(readFileSync(gboxFullPath, 'utf8')) as IGboxSpec;
        } else {
          return gbox;
        }
      })();

      /**
       * TODO: check the loaded in YAML against an auto-generated JSONSchema (from Typescript annotation)
       */
      await installOneGbox(knex, meta as object, gboxSpec);
    }),
  );
};
