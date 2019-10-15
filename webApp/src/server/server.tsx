import { createGenerateClassName, jssPreset } from '@material-ui/core/styles';
import { InMemoryCache } from 'apollo-cache-inmemory';
import { ApolloClient } from 'apollo-client';
import bodyParser from 'body-parser';
import chalk from 'chalk';
import cookieParser from 'cookie-parser';
import express from 'express';
import graphqlHTTP from 'express-graphql';
import expressJwt, { UnauthorizedError as Jwt401Error } from 'express-jwt';
import RateLimit from 'express-rate-limit';
import { ensureDir } from 'fs-extra';
import jwt from 'jsonwebtoken';
import * as Jss from 'jss';
import multer from 'multer';
import { resolve } from 'path';
import now from 'performance-now';
import { createPostGraphileSchema, mixed, withPostGraphileContext } from 'postgraphile';
import React from 'react';
import { getDataFromTree } from 'react-apollo';
import ReactDOMServer from 'react-dom/server';
import { SheetsRegistry } from 'react-jss/lib/jss';
import fetch from 'universal-fetch';
import uuidv4 from 'uuid/v4';

import App from '../common/App';
import configureStore from '../common/redux/configureStore';
// @ts-ignore
import assets from './assets.json';
import config from './config';
import createSandbox from './createSandbox';
import getPostgraphilePlugins from './getPostgraphilePlugins';
import getSingletonKnex from './getSingletonKnex';
import Html from './Html';
import startHttpsServer from './httpsServer';
import LocalLink from './LocalLink';
import passport from './passport';
import pgPool from './pgPool';
import queryProcessor from './queryProcesser';
// import models from './sequelize/models/index';

const upload = multer({
  storage: multer.diskStorage({
    destination: config.uploadPath,
    filename: (req, res, cb) => {
      cb(null, uuidv4());
    },
  }),
});
const app = express();

setImmediate(async () => {
  await ensureDir(config.dataPath);
  await ensureDir(config.downloadFilePath);
  await ensureDir(config.uploadPath);

  const knex = await getSingletonKnex();

  const postgraphilePlugins = await getPostgraphilePlugins(knex);
  const schema = await createPostGraphileSchema(config.databaseUrl, ['public'], {
    disableDefaultMutations: true,
    jwtPgTypeIdentifier: 'public.jwt_token',
    jwtSecret: config.auth.jwt.secret,
    prependPlugins: postgraphilePlugins,
    dynamicJson: true,
  });

  app.use(cookieParser());

  app.use(express.static(resolve(__dirname, 'public')));

  app.use(bodyParser.urlencoded({ extended: true }));
  app.use(bodyParser.json({ limit: '10gb' }));

  app.use('/sandbox', async (req, res) => {
    const sandboxId = await createSandbox(knex);
    // tslint:disable-next-line:no-console
    console.log(
      [' ii ', chalk.bold.red.inverse(` LOGIN `), ` The client has been assigned the sandbox id ${sandboxId}`].join(''),
    );
    const expiresIn = 60 * 60 * 24 * 180; // 180 days
    const token = jwt.sign(
      {
        role: 'logged_in_user',
        aud: ['postgraphile'],
        user_id: sandboxId,
      },
      config.auth.jwt.secret,
      { expiresIn },
    );
    res.cookie('id_token', token, { maxAge: 1000 * expiresIn, httpOnly: true });
    res.redirect('/');
  });

  app.get('/login', async (req, res) => {
    if (__DEV__) {
      res.redirect('/sandbox');
    } else {
      res.redirect('/sandbox');
    }
  });

  app.get('/logout', (req, res) => {
    res.clearCookie('id_token');
    res.redirect('/');
  });

  app.get('/login/debug', (req, res) => {
    const expiresIn = 60 * 60 * 24 * 180; // 180 days
    const token = jwt.sign(req.query, config.auth.jwt.secret, { expiresIn });
    res.cookie('id_token', token, { maxAge: 1000 * expiresIn, httpOnly: true });
    res.redirect('/');
  });

  app.use(
    async (req, res, next) => {
      if (req.cookies.id_token == null || req.cookies.id_token === '') {
        // tslint:disable-next-line:no-console
        console.log(
          [
            ' xx ',
            chalk.bold.red.inverse(` LOGIN `),
            ' The request does not have a "id_token" field in its cookies. Redirecting to /login ...',
          ].join(''),
        );
        res.redirect('/login');
        return;
      }

      next();
    },

    expressJwt({
      secret: config.auth.jwt.secret,
      credentialsRequired: false,
      getToken: (req) => req.cookies.id_token,
    }),

    async (req, res, next) => {
      if (req.user.user_id == null) {
        // tslint:disable-next-line:no-console
        console.log(
          [
            ' xx ',
            chalk.bold.red.inverse(` LOGIN `),
            ' The user does not have a user_id. Redirecting to /login ...',
          ].join(''),
        );
        res.redirect('/login');
        return;
      }

      const resCurrentUserAccount = await knex('user_account')
        .select(['id'])
        .where({ id: req.user.user_id });

      if (resCurrentUserAccount[0] == null) {
        // tslint:disable-next-line:no-console
        console.log(
          [
            ' xx ',
            chalk.bold.red.inverse(` LOGIN `),
            ' The user_id provided (',
            JSON.stringify(req.user.user_id),
            ') is not found in the database. Redirecting to /login ...',
          ].join(''),
        );
        res.redirect('/login');
        return;
      }

      next();
    },
  );

  // Error handler for express-jwt
  app.use((err, req, res, next) => {
    if (err instanceof Jwt401Error) {
      console.error('[express-jwt-error]', req.cookies.id_token);
      // `clearCookie`, otherwise user can't use web-app until cookie expires
      res.clearCookie('id_token');
    }
    next(err);
  });

  app.use(passport.initialize());

  if (__DEV__) {
    app.enable('trust proxy');
  }

  // app.get(
  //   '/login/facebook',
  //   passport.authenticate('facebook', {
  //     scope: ['email', 'user_location'],
  //     session: false,
  //   }),
  // );
  // app.get(
  //   '/login/facebook/return',
  //   passport.authenticate('facebook', {
  //     failureRedirect: '/login',
  //     session: false,
  //   }),
  //   (req, res) => {
  //     const expiresIn = 60 * 60 * 24 * 180; // 180 days
  //     const token = jwt.sign(req.user, config.auth.jwt.secret, { expiresIn });
  //     res.cookie('id_token', token, { maxAge: 1000 * expiresIn, httpOnly: true });
  //     res.redirect('/');
  //   },
  // );

  // const graphqlMiddleware = graphqlHTTP(({ user }) => ({
  //   schema,
  //   graphiql: __DEV__,
  //   context: { user },
  //   pretty: __DEV__,
  // }));

  // app.use('/graphql', graphqlMiddleware);

  // const getJWTToken = (req, res) => jwt.sign(req.user, config.auth.jwt.secret);
  //
  // const populateJWT = (req, res, next) => {
  //   if (req.headers.authorization == null) {
  //     req.headers.authorization = `bearer ${getJWTToken(req, res)}`;
  //   }
  //   next();
  // };

  app.post('/file-upload', upload.single('file'), (req, res) => {
    res.json({ fileId: req.file.filename });
  });

  // app.use('/graphql', populateJWT);
  // app.use('/graphiql', populateJWT);
  // app.use(
  //   postgraphile(config.databaseUrl, {
  //     graphiql: __DEV__,
  //     graphiqlRoute: '/graphiql',
  //     graphqlRoute: '/graphql',
  //
  //     // for createPostGraphQLSchema
  //     disableDefaultMutations: true,
  //     jwtPgTypeIdentifier: 'public.jwt_token',
  //     dynamicJson: true,
  //
  //     // for withPostGraphileContext
  //     jwtSecret: config.auth.jwt.secret,
  //     jwtAudiences: ['postgraphile'],
  //     pgDefaultRole: 'anonymous',
  //
  //     // for createPostGraphQLHttpRequestHandler
  //     showErrorStack: __DEV__,
  //     extendedErrors: __DEV__ ? ['hint', 'detail', 'errcode'] : [],
  //   }),
  // );

  app.use('/graphql', (req, res, next) => {
    withPostGraphileContext(
      {
        pgPool,
        jwtToken: req.cookies.id_token,
        jwtSecret: config.auth.jwt.secret,
        jwtRole: ['role'],
        jwtAudiences: ['postgraphile'],
        pgDefaultRole: 'anonymous',
      },
      // The reason why the following is a *callback* as opposed to the whole thing being a Promise is because,
      // after this callback is executed, there's one more query to make: "transaction commit".
      // That is, the callback is sandwiched between context-setting queries and transaction completion queries
      // This is why graphqlHTTP (which is async) has to be awaited.
      async (context: object) => {
        const start = now();

        await graphqlHTTP({
          schema,
          context: { ...context, user: req.user },
          graphiql: __DEV__,
          pretty: __DEV__,
        })(req, res);

        const end = now();
        // tslint:disable-next-line:no-console
        console.log(
          ` .. ` +
            `${chalk.bold.blue.inverse(` APOLLO QUERY `)} ` +
            `(${(end - start).toFixed(2)} ms) ` +
            `${req.body.operationName} `,
        );

        return null;
      },
    );
  });

  app.get('/download-data/:exportId', async (req, res) => {
    const exportItem = (await knex('export')
      .join('step', 'step.id', 'export.step_id')
      .join('project', 'project.id', 'step.project_id')
      .select(['export.id', 'export.extract_from', 'owner_id'])
      .where({ 'export.id': req.params.exportId }))[0];

    if (exportItem == null) {
      res.status(404).json({ error: 'File not found.' });
    }

    if (req.user.user_id !== exportItem.owner_id) {
      res.status(403).json({ error: 'You have no access to this file.' });
    }

    res
      .set({
        'Content-Disposition': `attachment; filename="${exportItem.extract_from}"`,
      })
      .sendFile(resolve(config.dataPath, exportItem.id));
  });

  app.post(
    '/query',
    new RateLimit({
      windowMs: 10 * 1000,
      max: 10,
      delayMs: 0,
      handler: (req, res, next) => {
        res.json({ errors: 'Too fast!!' });
      },
    }),
    await queryProcessor(),
  );

  app.get('*', async (req, res, next) => {
    // tslint:disable-next-line:no-console
    console.log(` <- a new request to ${req.path} (user = ${JSON.stringify(req.user)})`);
    try {
      // The schema has already been created and we are just retrieving the cache
      const apolloClient = new ApolloClient({
        ssrMode: true,
        cache: new InMemoryCache(),
        link: new LocalLink({
          schema,
          rootValue: {},
          contextValue: { jwtToken: req.cookies.id_token, user: req.user },
        }),
      });

      const initialState = {};

      const store = await configureStore(initialState, {
        apolloClient,
        fetch,
        req,
      });

      const location = store.getState().location;
      // At this point, the path thunks have all been executed,
      // all data related to the path have been fetched locally,
      // and they have resided in the store. (e.g. currentProject)

      if (location.kind === 'redirect') {
        res.redirect(302, location.pathname);
        // tslint:disable-next-line:no-console
        console.log(` -> redirected to ${location.pathname}`);
        return;
      }

      const sheetsRegistry = new SheetsRegistry();
      const jss = Jss.create(jssPreset() as any);
      const generateClassName = createGenerateClassName();

      const appProps = {
        userAgent: req.headers['user-agent'],
        store,
        generateClassName,
        sheetsRegistry,
        jss,
        apolloClient,
      };

      const rootComponent = <App {...appProps} />;
      await getDataFromTree(rootComponent);

      console.log('assets =', assets);

      const htmlProps = {
        title: store.getState().title,
        description:
          'Granatum is a graphical single-cell RNA-seq (scRNA-seq) analysis pipeline for genomics scientists.',
        preloads: [assets.vendor.js, assets.client.js],
        scripts: [assets.vendor.js, assets.client.js],
        app: {
          state: store.getState(),
          apollo: apolloClient.cache.extract(),
        },
        useLoadingScreen: true,
      };

      const html = ReactDOMServer.renderToStaticMarkup(<Html {...htmlProps} />);
      res.status(200);
      res.send(`<!doctype html>${html}`);
      // tslint:disable-next-line:no-console
      console.log(` -> sent 200`);
    } catch (err) {
      next(err);
    }
  });

  /**
   * This "afterware" ensures that even when a big error occurs at the server:
   *
   *   - The app itself won't die.
   *   - A graceful message is return to the user. (Debug messages, if in DEV mode)
   *
   */
  app.use(async (err, req, res, _next) => {
    console.error(err.stack);

    // TODO: make it look nicer
    const gracefulMessage = <pre>We are sorry. A server error just happend.</pre>;
    const debugMessage = (
      <div>
        <h1>Internal Error</h1>
        <pre>{err.stack}</pre>
      </div>
    );
    const html = ReactDOMServer.renderToStaticMarkup(
      <Html
        title="Internal Server Error"
        description={err.message}
        childrenString={ReactDOMServer.renderToString(__DEV__ ? debugMessage : gracefulMessage)}
      />,
    );
    res.status(err.status || 500);
    res.send(`<!doctype html>${html}`);
  });

  // Launch the server: this is for `npm run build`
  // -----------------------------------------------------------------------------
  if ((module as any).hot == null) {
    app.listen(config.port, '0.0.0.0', () => {
      console.info(`The server is running at http://localhost:${config.port}/`);
    });
  }

  // Set process title to make it easy to check running status of this app from other apps
  // Todo: process.title is limited in length to length of initial command, switch to a better mechanism
  process.title = 'granatumWebApp';

  // Try to do SSL even if only one var is set so a user gets a helpful error
  if (process.env.SSL_CERT || process.env.SSL_KEY) {
    startHttpsServer(app, process.env.SSL_CERT, process.env.SSL_KEY);
  }
});

// Hot Module Replacement
// -----------------------------------------------------------------------------
if ((module as any).hot) {
  (app as any).hot = (module as any).hot;
}

export default app;
