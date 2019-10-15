import chalk from 'chalk';
import { createLogger as reduxLogger } from 'redux-logger';
import { inspect } from 'util';

function inspectObject(object) {
  return inspect(object, {
    depth: 2,
    maxArrayLength: 4,
    colors: true,
    breakLength: Infinity,
  });
}

export default (process.env.BROWSER
  ? () => reduxLogger({ collapsed: true })
  : () => (_store) => (next) => (action) => {
      // tslint:disable-next-line:no-console
      console.log(` ** ${chalk.bold.inverse(` ${action.type} `)} ${inspectObject(action.payload)}`);
      return next(action);
    });
