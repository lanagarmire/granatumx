// @ts-ignore
import normalizeCss from 'normalize.css';

import React from 'react';
import serialize from 'serialize-javascript';

// eslint-disable-next-line no-unused-vars
// const GOOGLE_TRACK_SCRIPTS = [
//   config.analytics.googleTrackingId && (
//     <script
//       dangerouslySetInnerHTML={{
//         __html:
//           'window.ga=function(){ga.q.push(arguments)};ga.q=[];ga.l=+new Date;' +
//           `ga('create','${config.analytics.googleTrackingId}','auto');ga('send','pageview')`,
//       }}
//     />
//   ),
//   config.analytics.googleTrackingId && <script src="https://www.google-analytics.com/analytics.js" async defer />,
// ];

const Html: React.FunctionComponent<any> = ({
  title,
  description,
  preloads,
  scripts,
  app,
  childrenString,
  useLoadingScreen,
}) => (
  <html className="no-js" lang="en">
    <head>
      <meta charSet="utf-8" />
      <title>{title}</title>
      <meta name="description" content={description} />
      <meta httpEquiv="x-ua-compatible" content="ie=edge" />
      <style id="normalize" dangerouslySetInnerHTML={{ __html: normalizeCss }} />
      <meta name="viewport" content="width=device-width, initial-scale=1" />
      <link rel="apple-touch-icon" href="../../public/apple-touch-icon.png" />
      <link href="https://fonts.googleapis.com/css?family=Roboto+Mono|Roboto:300,400,500" rel="stylesheet" />
      <link href="https://fonts.googleapis.com/icon?family=Material+Icons" rel="stylesheet" />
      {preloads && preloads.map((preload) => <link key={preload} rel="preload" href={preload} as="script" />)}
    </head>
    <body>
      {childrenString && <div dangerouslySetInnerHTML={{ __html: childrenString }} />}
      <div id="app" />
      {useLoadingScreen && (
        <div
          style={{
            // display: 'none',
            zIndex: 99999,
            height: '100vh',
            width: '100vw',
            position: 'absolute',
            left: '0',
            top: '0',
            backgroundColor: '#fff',
          }}
          id="loading-blocker"
        >
          <div
            style={{
              position: 'fixed',
              top: '50%',
              left: '50%',
              transform: 'translate(-50%, -50%)',
              fontSize: '4em',
              fontFamily: 'sans',
            }}
          >
            Loading ...
          </div>
        </div>
      )}
      {app && <script dangerouslySetInnerHTML={{ __html: `window.App=${serialize(app)}` }} />}
      {scripts && scripts.map((script) => <script key={script} src={script} />)}
    </body>
  </html>
);

export default Html;
