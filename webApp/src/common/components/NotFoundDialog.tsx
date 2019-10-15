import { Button, Typography as T } from '@material-ui/core';
import React from 'react';
import { connect } from 'react-redux';

const NotFoundDialog = ({ prevLocation }) => (
  <div style={{ position: 'relative', width: '100vw', height: '100vh' }}>
    <div style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)' }}>
      <T variant="h2">404 Not Found</T>
      <T>Sorry ... we can&apos;t find what you are looking for.</T>

      <Button
        component={(props) => (
          // eslint-disable-next-line jsx-a11y/anchor-has-content
          <a href="/" {...props as any} />
        )}
      >
        Go to Homepage
      </Button>
      {prevLocation && (
        <Button
          component={(props) => (
            // eslint-disable-next-line jsx-a11y/anchor-has-content
            <a href={prevLocation} {...props as any} />
          )}
          color="primary"
        >
          Go back
        </Button>
      )}
    </div>
  </div>
);

export default connect((state: IReduxState) => ({
  prevLocation: state.location && state.location.prev && state.location.prev.pathname,
}))(NotFoundDialog);
