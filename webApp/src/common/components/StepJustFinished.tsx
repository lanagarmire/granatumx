import { Button, IconButton, Snackbar, withStyles } from '@material-ui/core';
import { Close } from '@material-ui/icons';
import gql from 'graphql-tag';
import React from 'react';
import { graphql } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, renderNothing, withProps } from 'recompose';

import { RFR_PROJECT_STEP, UPDATE_STEP_JUST_CHANGED_STATUS } from '../constants';
import { guardEmpty } from '../utils';
import actionCreators from '../redux/actionCreators';

const styles = {};

const StepJustFinished: React.FunctionComponent<any> = ({
  stepId,
  stepRank,
  stepTitle,
  rfrProjectStep,
  stepJustFinished,
  updateStepJustChangedStatus,
  classes,
}) => (
  <Snackbar
    open={stepJustFinished.open}
    anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
    autoHideDuration={6000}
    onClose={(e, reason) => {
      if (reason === 'clickaway') {
        return;
      }
      updateStepJustChangedStatus({ open: { $set: false } });
    }}
    ContentProps={{
      'aria-describedby': 'message-id',
    }}
    message={
      <span id="message-id">
        Step {stepRank + 1} <strong>{stepTitle}</strong> just finished.
      </span>
    }
    action={
      <>
        <Button
          color="secondary"
          size="small"
          onClick={() => {
            rfrProjectStep({ stepId });
            updateStepJustChangedStatus({ open: { $set: false } });
          }}
        >
          View results
        </Button>
        <IconButton
          key="close"
          aria-label="Close"
          color="inherit"
          className={classes.close}
          onClick={() => {
            updateStepJustChangedStatus({ open: { $set: false } });
          }}
        >
          <Close />
        </IconButton>
      </>
    }
  />
);

const enhance: any = compose(
  connect(
    (state: IReduxState) => ({
      stepJustFinished: state.app.dialog.stepJustFinished,
    }),
    {
      rfrProjectStep: actionCreators.rfrProjectStep,
      updateStepJustChangedStatus: actionCreators.updateStepJustChangedStatus,
    },
  ),
  withProps(({ stepJustFinished }) => ({
    stepId: stepJustFinished.stepId,
  })),
  branch(({ stepId }) => stepId == null, renderNothing),
  graphql(gql`
    query StepJustFinished($stepId: UUID!) {
      stepById(id: $stepId) {
        id
        rank
        gboxByGbox {
          id
          meta
        }
      }
    }
  `),
  branch(({ data }) => data.loading, renderNothing),
  guardEmpty('stepById'),
  withProps(({ data }) => ({ stepRank: data.stepById.rank, stepTitle: data.stepById.gboxByGbox.meta.title })),
  withStyles(styles),
);

export default enhance(StepJustFinished);
