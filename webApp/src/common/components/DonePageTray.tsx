import { Button, withStyles } from '@material-ui/core';
import React from 'react';
import { compose } from 'recompose';

const styles = (theme) => ({
  button: {
    margin: theme.spacing.unit,
    '&:first-child': {
      marginLeft: 0,
    },
    '&:last-child': {
      marginRight: 0,
    },
  },
});

const DonePageTray: React.FunctionComponent<any> = ({ gReset, gNextStep, classes }) => (
  <div>
    <Button variant="contained" color="primary" onClick={gReset} className={classes.button}>
      Reset
    </Button>
    <Button variant="contained" color="primary" onClick={gNextStep} className={classes.button}>
      Next step
    </Button>
  </div>
);

const enhance = compose(withStyles(styles));

export default enhance(DonePageTray);
