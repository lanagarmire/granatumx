import { ensureDirSync } from 'fs-extra';
import Knex from 'knex';
import config from './config';
import searchForTaskToRun from './searchForTaskToRun';
import clearRunningTasks from './clearRunningTasks';

// Set process title to make it easier to check taskRunner's running status from other apps
process.title = 'node taskRunner';

/**
 * TODO: Chores:
 *
 *   - remove all the transient files periodically,
 *   - remove all expired sandboxes
 *
 */

/**
 * TESTING
 * =======
 *
 * Preparing a simulated environment for testing a gbox:
 *
 *    - put the step data (`args.json`, `imports/`, `exports/` etc) in a folder (e.g. /a/b/example)
 *    - run the gbox with `env GRANATUM_SWD=/a/b/example docker run ...`
 *    - check errors in that directory
 *
 */

// TODO: needs to *not quit* when an error happens

setImmediate(async () => {
  ensureDirSync(config.uploadPathOnHost);
  ensureDirSync(config.dataPathOnHost);
  ensureDirSync(config.stepPathBaseOnHost);

  const knex = await Knex({
    client: 'pg',
    connection: process.env.DATABASE_URL || 'postgres://postgres:12qw@localhost:5433/granatum',
  });

  if (config.clearRunningTasks) {
    console.log('Try to clear any running tasks.');
    await clearRunningTasks(knex);
  }

  await searchForTaskToRun(knex);
});
