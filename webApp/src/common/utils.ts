import update from 'immutability-helper';
import { branch, compose, renderNothing, withHandlers, withState } from 'recompose';

export const range0 = (n: number) => {
  const output = [];
  for (let i = 0; i < n; i++) {
    output.push(i);
  }
  return output;
};

export const colsToRows = (cols: Record<string, any[]>): Array<Record<string, any>> => {
  const colnames = Object.keys(cols);
  if (colnames[0] == null) {
    return [];
  }
  return cols[colnames[0]].map((col, i) => {
    const row = {};
    colnames.forEach((n) => {
      row[n] = cols[n][i];
    });
    return row;
  });
};

export const trace = <T>(label: string, x: T): T => {
  // tslint:disable-next-line:no-console
  console.log(label, '=', x);
  return x;
};

export const defaultTo = (x) => (y) => (y == null ? x : y);

export const downloadURI = (uri, name = null) => {
  if (typeof document === 'object') {
    const link = document.createElement('a');
    if (name != null) {
      link.download = name;
    }
    link.href = uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  }
};

export const delay = (t) =>
  new Promise((pYes) => {
    setTimeout(() => {
      pYes();
    }, t);
  });

export const addId = (x) =>
  x && typeof x.map === 'function' ? x.map((y) => ({ _id: Math.random().toString(), ...y })) : x;

// noinspection CommaExpressionJS
export const guardEmpty = (dataName, propName = 'data') =>
  branch(
    (props) =>
      typeof props[propName][dataName] === 'undefined' &&
      (setTimeout(() => {
        if (__DEV__) {
          // tslint:disable-next-line:no-console
          console.error('guardEmpty');
        }
        props[propName].refetch();
      }),
      true),
    renderNothing,
  );

export const capitalize = (s) => s.slice(0, 1).toUpperCase() + s.slice(1);

export const withStateUpdater = (stateName, initialState) =>
  compose(
    withState(stateName, `set${capitalize(stateName)}`, initialState),
    withHandlers({
      [`update${capitalize(stateName)}`]: (props) => (updater) => {
        props[`set${capitalize(stateName)}`](update(props[stateName], updater));
      },
    }),
  );

export const queryBackend = async (fetch, payload) => {
  const response = await fetch('/query', {
    method: 'POST',
    body: JSON.stringify(payload),
    credentials: 'include',
    headers: new Headers({
      'Content-Type': 'application/json',
    }),
  });

  return response.json();
};
