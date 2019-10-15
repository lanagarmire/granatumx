import Knex from 'knex';

export default async (knex: Knex) => {
  const runningTasks: Array<{ id: string }> = await knex('step')
    .select(['id'])
    .where('status', 'running')
    .orWhere('status', 'interception_requested');

  if (runningTasks.length > 0) {
    console.log(`Found ${runningTasks.length} running tasks, resetting them ...`);
    await Promise.all(runningTasks.map((task) => knex.raw('select reset_step_recursively(?)', [task.id])));
    console.log(`Done.`);
  }
};
