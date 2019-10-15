import getSingletonKnex from './getSingletonKnex';

describe('test knex singleton', async () => {
  // note: jest seems to break if awaits are not inside the it tests
  let knex;
  it('should connect', async () => {
    knex = await getSingletonKnex();
    expect(knex).toBeTruthy();
  });

  it('should be singleton', async () => {
    const knex2 = await getSingletonKnex();
    expect(knex2).toBe(knex);
  });
});
