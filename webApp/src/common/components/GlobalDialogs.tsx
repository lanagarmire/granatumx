import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Tooltip,
  Typography as T,
  withStyles,
} from '@material-ui/core';
import { Error as ErrorIcon, Info as InfoIcon, Warning as WarningIcon } from '@material-ui/icons';
import React from 'react';
import ReactMarkdown from 'react-markdown';
import { connect } from 'react-redux';
import { compose, withProps } from 'recompose';

import { blue, orange, red } from '@material-ui/core/colors';
import { SvgIconProps } from '@material-ui/core/SvgIcon';
import actionCreators from '../redux/actionCreators';

const styles = {};

const GlobalDialogs = ({
  globalDialogs,
  dynamicDispatchAction,
  dismissGlobalDialog,
  typeToIcon,
}: {
  globalDialogs: IStateGlobalDialog[];
  dismissGlobalDialog;
  dynamicDispatchAction;
  typeToIcon: Record<string, { icon: React.ComponentClass<SvgIconProps>; color: string }>;
}) => (
  <>
    {globalDialogs.map((e) => (
      <Dialog key={e.id} open>
        {(e.type || e.title) && (
          <DialogTitle>
            {e.type
              ? React.createElement(typeToIcon[e.type].icon, {
                  style: {
                    marginRight: '0.2em',
                    marginBottom: '-0.2em',
                    color: typeToIcon[e.type].color,
                  },
                })
              : null}
            {e.title}
          </DialogTitle>
        )}
        <DialogContent>
          {typeof e.content === 'string' ? (
            <T component={'div' as any}>
              <ReactMarkdown source={e.content || ''} />
            </T>
          ) : (
            e.content
          )}
        </DialogContent>
        <DialogActions>
          {e.actions &&
            e.actions.map((a) => (
              <Tooltip key={a.label} title={a.description}>
                <Button
                  onClick={() => {
                    dismissGlobalDialog({ id: e.id });
                    dynamicDispatchAction({ action: a.action });
                  }}
                >
                  {a.label}
                </Button>
              </Tooltip>
            ))}
          <Button
            onClick={() => {
              dismissGlobalDialog({ id: e.id });
            }}
          >
            OK
          </Button>
        </DialogActions>
      </Dialog>
    ))}
  </>
);

const enhance = compose<any, any>(
  connect(
    (state: IReduxState) => ({
      globalDialogs: state.globalDialogs,
    }),
    {
      dismissGlobalDialog: actionCreators.dismissGlobalDialog,
      dynamicDispatchAction: actionCreators.dynamicDispatchAction,
    },
  ),
  withProps({
    typeToIcon: {
      error: {
        icon: ErrorIcon,
        color: red[600],
      },
      info: {
        icon: InfoIcon,
        color: blue[600],
      },
      warning: {
        icon: WarningIcon,
        color: orange[600],
      },
    },
  }),
  withStyles(styles),
);

export default enhance(GlobalDialogs);
