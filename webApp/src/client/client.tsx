/**
 * React Starter Kit (https://www.reactstarterkit.com/)
 *
 * Copyright © 2014-present Kriasoft, LLC. All rights reserved.
 * Copyright © 2018-present Xun Zhu. All rights reserved.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE.txt file in the root directory of this source tree.
 */

import { createGenerateClassName, jssPreset } from '@material-ui/core/styles';
// eslint-disable-next-line import/no-unresolved,import/extensions
import { InMemoryCache } from 'apollo-cache-inmemory';
import ApolloClient from 'apollo-client';
import { HttpLink } from 'apollo-link-http';
import { fromJS } from 'immutable';
import * as Jss from 'jss';
import React from 'react';
import ReactDOM from 'react-dom';
import { SheetsRegistry } from 'react-jss/lib/jss';
import fetch from 'universal-fetch';

import App from '../common/App';
import configureStore from '../common/redux/configureStore';
import pollingTick from './pollingTick';

(async () => {
  try {
    const apolloClient = new ApolloClient({
      cache: new InMemoryCache().restore((window as any).App.apollo),
      link: new HttpLink({
        uri: '/graphql',
        credentials: 'same-origin',
      }),
    });

    const container = document.getElementById('app');

    const jssStyles = document.getElementById('jss-server-side');
    if (jssStyles && jssStyles.parentNode) {
      jssStyles.parentNode.removeChild(jssStyles);
    }

    const sheetsRegistry = new SheetsRegistry();
    const jss = Jss.create(jssPreset() as any);
    const generateClassName = createGenerateClassName();

    // TODO: make the entire state ImmutableJS
    (window as any).App.state.app.componentStates = fromJS((window as any).App.state.app.componentStates);

    const store = await configureStore((window as any).App.state, {
      apolloClient,
      fetch,
    });

    store.subscribe(() => {
      document.title = store.getState().title;
    });

    const appProps = {
      userAgent: navigator.userAgent,
      store,
      jss,
      sheetsRegistry,
      generateClassName,
      apolloClient,
    };

    // The polling tick
    pollingTick(store, apolloClient);

    ReactDOM.render(<App {...appProps} />, container);

    setTimeout(() => {
      const loadingBlocker = document.getElementById('loading-blocker');
      if (loadingBlocker && loadingBlocker.parentNode) {
        loadingBlocker.parentNode.removeChild(loadingBlocker);
      }
    });
  } catch (error) {
    if (__DEV__) {
      throw error;
    }

    console.error(error);
  }
})();
