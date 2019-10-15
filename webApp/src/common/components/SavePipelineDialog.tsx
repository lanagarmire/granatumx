import { Button, Dialog, DialogActions, DialogContent, DialogTitle, TextField } from '@material-ui/core';
import gql from 'graphql-tag';
import _ from 'lodash';
import React from 'react';
import { graphql } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, renderNothing, withProps } from 'recompose';
import { SAVE_CURRENT_PIPELINE_AS_RECIPE } from '../constants';

import { guardEmpty, withStateUpdater } from '../utils';
import Pipeline from './Pipeline';
import actionCreators from '../redux/actionCreators';

const SavePipelineDialog: React.FunctionComponent<any> = ({
  title,
  setTitle,
  titleError,
  setTitleError,
  subtitle,
  setSubtitle,
  subtitleError,
  saveCurrentPipelineAsRecipe,
  open,
  onClose,
  steps,
  ...props
}) => (
  <Dialog open={open} onClose={onClose}>
    <DialogTitle>Save pipeline as a recipe</DialogTitle>
    <DialogContent>
      <TextField
        fullWidth
        error={!!titleError}
        helperText={titleError}
        label="Recipe title"
        InputLabelProps={{
          shrink: true,
        }}
        value={title}
        onChange={(e) => {
          if (e.target.value.length === 0) {
            setTitleError('Empty title not allowed.');
          } else {
            setTitleError(null);
          }
          setTitle(e.target.value);
        }}
      />
    </DialogContent>
    <DialogContent>
      <TextField
        fullWidth
        error={!!subtitleError}
        helperText={subtitleError}
        label="Recipe subtitle"
        InputLabelProps={{
          shrink: true,
        }}
        value={subtitle}
        onChange={(e) => {
          setSubtitle(e.target.value);
        }}
      />
    </DialogContent>
    <DialogContent>
      <Pipeline steps={steps} />
    </DialogContent>
    <DialogActions>
      <Button
        onClick={() => {
          onClose();
        }}
      >
        Cancel
      </Button>
      <Button
        color="primary"
        onClick={() => {
          saveCurrentPipelineAsRecipe({ meta: { title, subtitle } });
          onClose();
        }}
      >
        Save
      </Button>
    </DialogActions>
  </Dialog>
);

const enhance: any = compose(
  connect(
    (state: IReduxState) => ({
      // currentProjectRank: state.location.payload.projectRank,
      currentProjectId: state.app.currentProjectId,
    }),
    {
      saveCurrentPipelineAsRecipe: actionCreators.saveCurrentPipelineAsRecipe,
    },
  ),
  branch((props: any) => props.currentProjectId == null, renderNothing),
  graphql(
    gql`
      query SavePipelineDialog($currentProjectId: UUID!) {
        projectById(id: $currentProjectId) {
          id
          stepsByProjectId {
            nodes {
              id
              rank
              status
              errors
              gboxByGbox {
                id
                meta
              }
            }
          }
        }
      }
    `,
  ),
  branch((props: any) => props.data.loading, renderNothing),
  guardEmpty('projectById'),
  withProps(({ data }) => ({
    steps: _(data.projectById.stepsByProjectId.nodes)
      .orderBy('rank')
      .map((x) => ({ id: x.id, title: x.gboxByGbox.meta.title }))
      .value(),
  })),
  branch(({ steps }) => steps == null, renderNothing),
  withStateUpdater('title', 'Recipe generated on ' + new Date().toUTCString()),
  withStateUpdater('titleError', null),
  withStateUpdater('subtitle', 'Saved from a custom pipeline'),
  withStateUpdater('subtitleError', null),
);

export default enhance(SavePipelineDialog);
