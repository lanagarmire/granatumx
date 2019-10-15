import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  Typography as T,
  withStyles,
  withWidth,
} from '@material-ui/core';
import cn from 'classnames';
import React from 'react';
import { compose } from 'recompose';
import { withStateUpdater } from '../utils';
import ReactMarkdown = require('react-markdown');

const styles: any = {
  root: {
    cursor: 'zoom-in',
  },
  dialogContent: {
    textAlign: 'center',
  },
  pngZoomDialogPaper: {
    width: '100%',
    height: '100%',
    maxWidth: '100%',
  },
  zoomedCaption: {
    flex: 1,
    padding: '0 16px',
  },
};

const ZoomableImage: React.FunctionComponent<any> = ({
  zoomDialog,
  setZoomDialog,
  classes,
  width,
  height,
  data,
  description,
  className,
}) => (
  <>
    <div
      onClick={() => {
        setZoomDialog(true);
      }}
      className={cn(classes.root, className)}
    >
      <img width={width} height={height} alt="plot" src={`data:image/png;base64,${data}`} />
      <T variant="caption">
        <ReactMarkdown source={description} />
      </T>
    </div>
    <Dialog
      classes={{ paper: classes.pngZoomDialogPaper }}
      open={zoomDialog}
      onClose={() => {
        setZoomDialog(false);
      }}
    >
      <DialogContent className={classes.dialogContent}>
        <img alt="plot" src={`data:image/png;base64,${data}`} />
      </DialogContent>
      <DialogActions>
        <div className={classes.zoomedCaption}>
          <T>
            <ReactMarkdown source={description} />
          </T>
        </div>
        <Button
          onClick={() => {
            setZoomDialog(false);
          }}
        >
          Close
        </Button>
      </DialogActions>
    </Dialog>
  </>
);

const enhance: any = compose(
  withStateUpdater('zoomDialog', false),
  withStyles(styles),
);

export default enhance(ZoomableImage);
