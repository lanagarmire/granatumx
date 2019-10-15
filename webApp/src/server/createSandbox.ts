export default async (knex) => {
  const userAccounts = await knex('user_account')
    .insert({ email: 'sandbox' })
    .returning('*');

  if (userAccounts[0] == null) {
    throw new Error('Sandbox creation was not successful.');
  }

  const sandboxId = userAccounts[0].id;

  await knex('user_profile').insert({ user_id: sandboxId, display_name: 'Sandbox' });

  return sandboxId;
};
