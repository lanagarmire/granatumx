import Knex from 'knex';
import { resolve } from 'path';
import config from './config';
// @ts-ignore
import stringArgv from 'string-argv';

export const getExeCmdArgv = async (
  knex: Knex,
  gboxId: string,
  taskRelatedFilesOnSWD: Array<{ path: string; opt?: string }>,
  swdOnHost: string,
) => {
  const gbox = (await knex('gbox')
    .select(['meta', 'endpoints'])
    .where('id', gboxId))[0];

  const endpointSpec = gbox.endpoints.backend;

  if (typeof endpointSpec === 'string') {
    return stringArgv(endpointSpec);
  }

  if (typeof endpointSpec === 'object') {
    if (endpointSpec.type === 'docker') {
      const cmd = [
        ...`docker run -i --rm -e GRANATUM_SWD=/data`.split(' '),
        ...taskRelatedFilesOnSWD
          .map((f) => ['-v', `${resolve(swdOnHost, f.path)}:${resolve(config.swdOnDocker, f.path)}:${f.opt}`])
          .reduce((a, x) => a.concat(x)),
        endpointSpec.image,
        ...stringArgv(endpointSpec.cmd),
      ];

      console.log('(for docker) cmd =', cmd);

      return cmd;
    }
  }

  throw new Error(`The spec ${JSON.stringify(endpointSpec)} is not recognized.`);
};
