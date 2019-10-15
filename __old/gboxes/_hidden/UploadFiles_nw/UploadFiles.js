import React from 'react';
import Dropzone from 'react-dropzone';
import {
  Grid,
  Typography as T,
  IconButton,
  TextField,
  Button,
  withStyles,
  Tooltip,
  Input,
  InputLabel,
  InputAdornment,
  FormControl,
  FormHelperText,
} from 'material-ui';
import { AttachFile, Replay, ChevronRight } from '@material-ui/icons';
import { compose, lifecycle, withHandlers, withProps, withState, withStateHandlers } from 'recompose';
import update from 'immutability-helper';
import jp from 'jsonpath';

const styles = {
  actionButtons: {
    position: 'fixed',
    bottom: 30,
    right: 350,
    '& > *': {
      margin: 10,
    },
  },
  pointer: {
    cursor: 'pointer',
  },
};

const initialState = {
  expressionMatrix: null,
  sampleMeta: null,
  assayName: 'Counts',
};

const UploadFiles = ({
  classes,
  onDropExpressionMatrix,
  onDropSampleMeta,
  onClickNext,
  gState = initialState,
  gSetState,
}) => (
  <div>
    <T variant="display4">Upload Files</T>
    <Grid container spacing={24}>
      <Grid item xs={12}>
        <Dropzone style={{}} onDrop={onDropExpressionMatrix} multiple={false}>
          <FormControl fullWidth>
            <InputLabel htmlFor="expressionMatrixFile">Expression Matrix</InputLabel>
            <Input
              classes={{ input: classes.pointer }}
              id="expressionMatrixFile"
              endAdornment={
                <InputAdornment>
                  <IconButton>
                    <AttachFile />
                  </IconButton>
                </InputAdornment>
              }
              value={gState.expressionMatrix ? gState.expressionMatrix.name : ''}
            />
          </FormControl>
        </Dropzone>
      </Grid>
      <Grid item xs={12}>
        <Dropzone style={{}} onDrop={onDropSampleMeta} multiple={false}>
          <FormControl fullWidth>
            <InputLabel htmlFor="sampleMetaFile">Sample Metadata</InputLabel>
            <Input
              classes={{ input: classes.pointer }}
              id="sampleMetaFile"
              endAdornment={
                <InputAdornment>
                  <IconButton>
                    <AttachFile />
                  </IconButton>
                </InputAdornment>
              }
              value={gState.sampleMeta ? gState.sampleMeta.name : ''}
            />
          </FormControl>
        </Dropzone>
      </Grid>
    </Grid>
    <div className={classes.actionButtons}>
      <Tooltip title="Reset the step">
        <Button fab onClick={() => gSetState(initialState)}>
          <Replay />
        </Button>
      </Tooltip>
      <Tooltip title="Continue to the next step">
        <Button fab color="primary" onClick={onClickNext}>
          <ChevronRight />
        </Button>
      </Tooltip>
    </div>
  </div>
);

const enhance = compose(
  withHandlers({
    onDropExpressionMatrix: (props) => (accepted) => {
      props.gSetState({ ...props.gState, expressionMatrix: accepted[0] });
    },
    onDropSampleMeta: (props) => (accepted) => {
      props.gSetState({ ...props.gState, sampleMeta: accepted[0] });
    },
    onClickNext: (props) => () => {
      console.log(props.gState);
    },
  }),
  withStyles(styles),
);

export const UPLOAD_FILES = 'UPLOAD_FILES';
export default enhance(UploadFiles);
