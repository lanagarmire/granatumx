/**
 * Because we are using "redux-first-router" pattern, the router configuration is
 * part of the store configuration.
 *
 */

import { applyMiddleware, combineReducers, compose, createStore } from 'redux';
import { connectRoutes } from 'redux-first-router';

import ApolloClient from 'apollo-client/ApolloClient';
import { Request } from 'express';
import createLogger from './logger';
import componentStates from './reducers/app/componentStates';
import currentProjectId from './reducers/app/currentProjectId';
import currentStepId from './reducers/app/currentStepId';
import addStep from './reducers/app/dialog/addStep';
import notFound from './reducers/app/dialog/notFound';
import stepJustFinished from './reducers/app/dialog/stepJustFinished';
import welcome from './reducers/app/dialog/welcome';
import globalDialogs from './reducers/globalDialogs';
import title from './reducers/title';
import routesMap from './routesMap';

export default async function configureStore(
  initialState,
  { apolloClient, fetch, req }: { apolloClient?: ApolloClient<any>; fetch?: any; req?: Request },
) {
  const rfr = connectRoutes(routesMap, {
    title: () => null,
    extra: { apolloClient, fetch },
    initialEntries: process.env.BROWSER ? undefined : req.path,
  });

  const middleware = [rfr.middleware];

  let enhancer;
  if (__DEV__) {
    middleware.push(createLogger());

    // https://github.com/zalmoxisus/redux-devtools-extension#redux-devtools-extension
    let devToolsExtension = (f) => f;
    if (process.env.BROWSER && (window as any).devToolsExtension) {
      devToolsExtension = (window as any).devToolsExtension();
    }

    enhancer = compose(
      applyMiddleware(...middleware),
      devToolsExtension,
    );
  } else {
    enhancer = applyMiddleware(...middleware);
  }

  enhancer = compose(
    rfr.enhancer,
    enhancer,
  );

  const rootReducer = combineReducers({
    location: rfr.reducer,
    app: combineReducers({
      componentStates,
      currentStepId,
      currentProjectId,
      dialog: combineReducers({
        welcome,
        notFound,
        addStep,
        stepJustFinished,
      }),
    }),
    title,
    globalDialogs,
  });

  // See https://github.com/rackt/redux/releases/tag/v3.1.0
  const store = createStore(rootReducer, initialState, enhancer);

  await rfr.thunk(store);

  return store;
}
