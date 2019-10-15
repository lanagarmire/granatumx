import { Button, LinearProgress, Typography as T, withStyles } from '@material-ui/core';
import React from 'react';
import { compose } from 'redux';
import escapeCarriage from 'escape-carriage';

const styles: any = {
  root: {
    position: 'fixed',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
  },
  actionTray: {
    margin: 16,
  },
  messageBox: {
    width: 750,
    maxHeight: 500,
    overflowX: 'auto',
    overflowY: 'auto',
  },
  message: {
    color: 'rgba(0, 0, 0, 0.5)',
    fontFamily: 'monospace',
  },
  messageHeader: {
    fontWeight: 'bold',
    color: 'rgba(0, 0, 0, 0.5)',
  },
};

const LoadingScreen = ({ label, message, actions, classes }) => (
  <div className={classes.root}>
    <div>
      <T variant="h2">{label}</T>
      <LinearProgress />
    </div>
    {message != null && (
      <div className={classes.messageBox}>
        {message.stdout == null || message.stdout.length === 0 ? null : (
          <>
            <T className={classes.messageHeader}>Stdout</T>
            <T component={'pre' as any} className={classes.message}>
              {escapeCarriage(message.stdout)}
            </T>
          </>
        )}
        {message.stderr == null || message.stderr.length === 0 ? null : (
          <>
            <T className={classes.messageHeader}>Stderr</T>
            <T component={'pre' as any} className={classes.message}>
              {escapeCarriage(message.stderr)}
            </T>
          </>
        )}
      </div>
    )}
    {actions && actions.length > 0 && (
      <div className={classes.actionTray}>
        {actions.map((x) => (
          <Button key={x.label} onClick={x.action}>
            {x.label}
          </Button>
        ))}
      </div>
    )}
  </div>
);

const enhance: any = compose(withStyles(styles));

export default enhance(LoadingScreen);
