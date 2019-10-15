import { ApolloLink, Observable } from 'apollo-link';
import chalk from 'chalk';
import { execute, ExecutionResult } from 'graphql';
import now from 'performance-now';
import { withPostGraphileContext } from 'postgraphile';

import config from './config';
import pgPool from './pgPool';

export default class LocalLink extends ApolloLink {
  private schema: any;
  private contextValue: any;
  private rootValue: any;

  constructor({ schema, rootValue, contextValue }) {
    super();

    this.schema = schema;
    this.contextValue = contextValue;
    this.rootValue = rootValue;
  }

  public request(operation) {
    const { schema, rootValue, contextValue } = this;
    const { query, variables, operationName } = operation;

    return new Observable((observer) => {
      let canceled = false;

      withPostGraphileContext(
        {
          pgPool,
          jwtToken: contextValue.jwtToken,
          jwtSecret: config.auth.jwt.secret,
          jwtRole: ['role'],
          jwtAudiences: ['postgraphile'],
          pgDefaultRole: 'anonymous',
        },
        async (context: object) => {
          const start = now();

          const result = await execute(
            schema,
            query,
            rootValue,
            { ...context, ...contextValue },
            variables,
            operationName,
          );

          const end = now();
          // tslint:disable-next-line:no-console
          console.log(
            ` .. ` +
              `${chalk.bold.blue.inverse(` APOLLO QUERY `)} ` +
              `(${(end - start).toFixed(2)} ms, ${JSON.stringify(result).length} bytes) ` +
              `${query.definitions.map((x) => x.name.value).join(', ')} `,
          );
          return result as ExecutionResult;
        },
      )
        .then((result) => {
          if (canceled) {
            return null;
          }
          // we have data and can send it to back up the link chain
          observer.next(result);
          observer.complete();

          return result;
        })
        .catch((err) => {
          if (canceled) {
            return;
          }

          observer.error(err);
        });

      return () => {
        canceled = true;
      };
    });
  }
}
