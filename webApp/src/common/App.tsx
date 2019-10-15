import { MuiThemeProvider } from '@material-ui/core';
import React from 'react';
import { ApolloProvider } from 'react-apollo';
import JssProvider from 'react-jss/lib/JssProvider';
import { Provider } from 'react-redux';

import { mainTheme } from './themes';

import LayoutChooser from './components/LayoutChooser';

const App = ({ store, sheetsRegistry, generateClassName, jss, apolloClient }) => (
  <Provider store={store}>
    <ApolloProvider client={apolloClient}>
      <JssProvider registry={sheetsRegistry} generateClassName={generateClassName} jss={jss}>
        <MuiThemeProvider theme={mainTheme} sheetsManager={new Map()}>
          <LayoutChooser />
        </MuiThemeProvider>
      </JssProvider>
    </ApolloProvider>
  </Provider>
);

export default App;
