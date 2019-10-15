import { Dialog, DialogContent, DialogTitle, TextField, withStyles, withWidth } from '@material-ui/core';
import React from 'react';
import { connect } from 'react-redux';
import { compose } from 'recompose';

import actionCreators from '../redux/actionCreators';
import { withStateUpdater } from '../utils';
import AddRecipe from './AddRecipe';
import AddStep from './AddStep';

const styles = {
  dialogPaper: {
    maxWidth: '100%',
    width: 800,
    height: '100%',
  },
};

const AddStepDialog: React.FunctionComponent<any> = ({
  updateAddStepDialog,
  addStepDialog,
  classes,
  filterText,
  setFilterText,
  width,
}) => (
  <Dialog
    onClose={() => {
      updateAddStepDialog({ open: { $set: false } });
    }}
    open={addStepDialog.open}
    classes={{ paper: classes.dialogPaper }}
    fullScreen={width === 'xs' || width === 'sm'}
  >
    <DialogContent style={{ flex: 'none' }}>
      <TextField
        fullWidth
        placeholder="Type to search ..."
        inputRef={(el) => {
          if (el != null) {
            el.focus();
          }
        }}
        onChange={(e) => {
          setFilterText(e.target.value);
        }}
        value={filterText}
      />
    </DialogContent>
    <div style={{ display: 'flex', flexDirection: 'row' }}>
      {addStepDialog.showRecipes && (
        <div style={{ flex: 1 }}>
          <AddRecipe
            projectId={addStepDialog.projectId}
            atStepRank={addStepDialog.atStepRank}
            filterText={filterText}
            onClose={() => {
              updateAddStepDialog({ open: { $set: false } });
            }}
          />
        </div>
      )}
      {addStepDialog.showSteps && (
        <div style={{ flex: 1 }}>
          <AddStep
            projectId={addStepDialog.projectId}
            atStepRank={addStepDialog.atStepRank}
            filterText={filterText}
            onClose={() => {
              updateAddStepDialog({ open: { $set: false } });
            }}
          />
        </div>
      )}
    </div>
  </Dialog>
);

const enhance = compose(
  connect(
    (state: IReduxState) => ({ addStepDialog: state.app.dialog.addStep }),
    {
      updateAddStepDialog: actionCreators.updateAddStepDialog,
    },
  ),
  withStateUpdater('filterText', ''),
  withWidth(),
  withStyles(styles),
);

export default enhance(AddStepDialog);
