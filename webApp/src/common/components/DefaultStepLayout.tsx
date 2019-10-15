import {
  Button,
  CircularProgress,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Typography as T,
  withStyles,
} from '@material-ui/core';
import { StyleRules } from '@material-ui/core/styles';
import React from 'react';
import ReactMarkdown from 'react-markdown';
import { compose, withState } from 'recompose';
import ZScatter from './ZScatter';

const styles: StyleRules<string> = {
  title: {
    marginBottom: '8px',
  },
  description: {
    marginBottom: '16px',
  },
  sectionTitle: {
    fontSize: '1.5em',
  },
  paper: {
    marginTop: '16px',
    padding: '16px',
  },
  form: {
    marginBottom: '32px',
    // '& > *:first-child': {
    //   marginTop: '0',
    // },
  },
  page: {
    margin: '16px 0 32px 0',
  },
  buttonRow: {
    margin: 0,
    marginLeft: 'auto',
    padding: 4,
    display: 'flex',
    flexDirection: 'row',
    backgroundColor: 'rgba(250, 250, 250, 0.7)',
    justifyContent: 'flex-end',
  },
  button: {
    margin: '4px',
    position: 'relative',
  },
  errorDialog: {
    maxWidth: 'none',
  },
};

const DefaultStepLayout = ({
  gErrors,
  gClearErrors,
  submitStep,
  disableSubmit,
  meta,
  children,

  showErrors,
  setShowErrors,
  classes,
}) => (
  <div className={classes.page}>
    <T variant="h2" gutterBottom>
      {meta.title}
    </T>
    <div className={classes.description}>
      <T gutterBottom component={'div' as any}>
        <ReactMarkdown source={meta.description} />
      </T>
      {meta.maintainer &&
        meta.maintainer.name && (
          <T variant="caption" align="right">
            Maintainer: {meta.maintainer.name}
          </T>
        )}
      {meta.maintainer &&
        meta.maintainer.email && (
          <T variant="caption" align="right">
            (<a href={`mailto:${meta.maintainer.email}`}>{meta.maintainer.email}</a>)
          </T>
        )}
    </div>
    {children}
    <div className={classes.buttonRow}>
      <Button
        className={classes.button}
        disabled={disableSubmit}
        variant="contained"
        color="primary"
        onClick={submitStep}
      >
        Submit
        {disableSubmit && (
          <CircularProgress
            size={24}
            style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)' }}
          />
        )}
      </Button>
    </div>
    <Dialog classes={{ paper: classes.errorDialog }} open={gErrors != null && showErrors}>
      <DialogTitle>The following error(s) occurred during the last execution</DialogTitle>
      {meta.maintainer &&
        meta.maintainer.name &&
        meta.maintainer.email && (
          <DialogContent>
            <T>
              If the problem persists, please contact the Gbox maintainer: {meta.maintainer.name} (
              <a href={`mailto:${meta.maintainer.email}`}>{meta.maintainer.email}</a>)
            </T>
          </DialogContent>
        )}
      <DialogContent>
        {typeof gErrors === 'string' ? (
          <T>
            <code>{gErrors}</code>
          </T>
        ) : typeof gErrors === 'object' && gErrors != null && typeof gErrors.map === 'function' ? (
          <Table>
            <TableHead>
              <TableRow>
                <TableCell>Source</TableCell>
                <TableCell>Message</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {gErrors.map((e, i) => (
                // eslint-disable-next-line react/no-array-index-key
                <TableRow key={i}>
                  <TableCell>{e.source}</TableCell>
                  <TableCell>
                    <T variant="caption">
                      <code>
                        <pre>{e.message}</pre>
                      </code>
                    </T>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        ) : (
          <T>
            <code>{JSON.stringify(gErrors)}</code>
          </T>
        )}
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => {
            gClearErrors();
            setShowErrors(false);
          }}
        >
          OK
        </Button>
      </DialogActions>
    </Dialog>
  </div>
);

// TODO: trying to use globalError instead of this ad hoc error reporting

const enhance = compose(
  withState('showErrors', 'setShowErrors', true),
  withStyles(styles),
);

export default enhance(DefaultStepLayout);
