import { Avatar, Button, Fab, ListItem, ListItemText, Tooltip, withStyles } from '@material-ui/core';
import { blue, common, green, purple, red, yellow } from '@material-ui/core/colors';
import { CSSProperties } from '@material-ui/core/styles/withStyles';
import { Add, Clear, Grain, Replay, UnfoldMore } from '@material-ui/icons';
import cn from 'classnames';
import React, { FunctionComponent } from 'react';
import { connect } from 'react-redux';
import { compose, withProps } from 'recompose';

import { REMOVE_STEP, RESET_STEP, RFR_NO_STEPS, UPDATE_ADD_STEP_DIALOG } from '../constants';
import actionCreators from '../redux/actionCreators';
import { withStateUpdater } from '../utils';

const styles: Record<string, CSSProperties> = {
  '@keyframes runningKeyframes': {
    '0%': { backgroundColor: '#fff' },
    '50%': { backgroundColor: green[100] },
    '100%': { backgroundColor: '#fff' },
  },
  stepNumIconLabel: {
    display: 'block',
    position: 'absolute',
    left: '50%',
    top: '50%',
    transform: 'translate(-50%,-50%)',
    transition: 'opacity 0.25s',
    opacity: 1,
  },
  stepNumIconLabelHover: {
    opacity: 0,
  },
  stepNumIconHandle: {
    color: '#777',
    display: 'block',
    position: 'absolute',
    left: '50%',
    top: '50%',
    transform: 'translate(-50%,-50%)',
    transition: 'opacity 0.25s',
    opacity: 0,
  },
  stepNumIconHandleHover: {
    color: '#fff',
    opacity: 1,
  },
  stepNumIcon: {
    margin: '16px 0 16px 16px',
    width: 30,
    height: 30,
    fontSize: '1rem',
    color: '#333',
    boxShadow: '0px 0px 2px 2px #eee',
    backgroundColor: '#fff',
  },
  stepNumIconHover: {
    backgroundColor: yellow[900],
  },
  stepNumIconCurrent: {
    // boxShadow: 'none',
  },
  itemContainer: {
    position: 'relative',
    width: '100%',
  },
  listItemButton: {
    paddingLeft: 0,
    paddingRight: 0,
  },
  listItem: {
    padding: 0,
    display: 'flex',
    flexDirection: 'row',
    backgroundColor: common.white,
  },
  listItemText: {
    padding: 16,
    flex: '1 1 0%',
  },
  listItemTextPrimary: {
    fontSize: '0.9rem',
  },
  listItemTextSecondary: {
    fontSize: '0.8rem',
  },
  done: {
    backgroundColor: green[50],
    '&:hover': {
      background: green[100],
    },
    '&:focus': {
      background: green[100],
    },
  },
  failed: {
    backgroundColor: red[100],
    '&:hover': {
      background: red[200],
    },
    '&:focus': {
      background: red[200],
    },
  },
  running: {
    animation: 'runningKeyframes 1s infinite',
  },
  current: {
    backgroundColor: purple[900],
    '&:hover': {
      background: purple[800],
    },
    '&:focus': {
      background: purple[800],
    },
  },
  currentText: {
    color: common.white,
  },
  dragging: {
    boxShadow: [
      '0px 11px 15px -7px rgba(0, 0, 0, 0.2)',
      '0px 24px 38px 3px rgba(0, 0, 0, 0.14)',
      '0px 9px 46px 8px rgba(0, 0, 0, 0.12)',
    ].join(),
  },
  functionButtonTray: {
    margin: '-4px',
    padding: '16px 16px 16px 0',
    position: 'absolute',
    zIndex: 1250, // above drawers, below modals
    bottom: -32,
    right: 0,
    display: 'flex',
    flexDirection: 'row',
    visibility: 'hidden',
  },
  functionButtonTrayHover: {
    visibility: 'visible',
  },
  functionButton: {
    width: 36,
    height: 36,
    margin: 4,
  },
  functionButtonIcon: {
    fontSize: '1rem',
    margin: 4,
  },
};

const StepListItem: FunctionComponent<any> = ({
  step,
  dragHandleProps,
  rfrProjectStep,
  currentProjectId,
  currentStepId,
  onScrollToStep,
  /**/
  isHovering,
  setIsHovering,
  functionButtons,
  classes,
}) => (
  <div
    onMouseOver={() => {
      setIsHovering(true);
    }}
    onMouseOut={() => {
      setIsHovering(false);
    }}
    className={classes.itemContainer}
  >
    <ListItem
      divider
      buttonRef={(el) => {
        if (step.id === currentStepId) {
          onScrollToStep({ stepElement: el });
        }
      }}
      button
      onClick={() => rfrProjectStep({ stepId: step.id })}
      classes={{
        button: classes.listItemButton,
        root: cn(
          classes.listItem,
          step.id === currentStepId
            ? classes.current
            : step.status === 'INITIATED' || step.status === 'RUNNING'
            ? classes.running
            : step.status === 'DONE'
            ? classes.done
            : step.errors != null
            ? classes.failed
            : null,
        ),
      }}
    >
      <div style={{ flex: 'none' }} {...dragHandleProps}>
        <Avatar
          className={cn(
            classes.stepNumIcon,
            {
              [classes.stepNumIconCurrent]: step.id === currentStepId,
            },
            isHovering && classes.stepNumIconHover,
          )}
        >
          <div className={cn(classes.stepNumIconLabel, isHovering && classes.stepNumIconLabelHover)}>
            {step.rank + 1}
          </div>
          <Tooltip title="Drag to move the step in the pipeline. Drag outside to remove the step">
            <UnfoldMore className={cn(classes.stepNumIconHandle, isHovering && classes.stepNumIconHandleHover)} />
          </Tooltip>
        </Avatar>
      </div>
      <ListItemText
        primary={`${step.gboxByGbox.meta.title}${step.status === 'INITIATED' || step.status === 'RUNNING' ? ' …' : ''}${
          step.status === 'DONE' ? ' ✓' : ''
        }`}
        secondary={step.gboxByGbox.meta.subtitle}
        classes={{
          root: classes.listItemText,
          primary: cn(classes.listItemTextPrimary, step.id === currentStepId && classes.currentText),
          secondary: cn(classes.listItemTextSecondary, step.id === currentStepId && classes.currentText),
        }}
      />
    </ListItem>
    <div className={cn(classes.functionButtonTray, isHovering && classes.functionButtonTrayHover)}>
      {functionButtons.map((f) => (
        <Tooltip key={f.label} title={f.label}>
          <Fab
            color={f.color}
            onClick={(e) => {
              // e.stopPropagation();
              f.onClick(e);
            }}
            className={classes.functionButton}
          >
            {React.createElement(f.icon, { className: classes.functionButtonIcon })}
          </Fab>
        </Tooltip>
      ))}
    </div>
  </div>
);

const enhance = compose(
  connect(
    null,
    {
      updateAddStepDialog: actionCreators.updateAddStepDialog,
      removeStep: actionCreators.removeStep,
      resetStep: actionCreators.resetStep,
    },
  ),
  withProps(({ resetStep, currentStepId, currentProjectId, removeStep, updateAddStepDialog, step }) => ({
    functionButtons: [
      ...(step.status !== 'DONE'
        ? []
        : [
            {
              label: 'Reset',
              color: 'default',
              icon: Replay,
              onClick: async () => {
                resetStep({ stepId: step.id });
              },
            },
          ]),
      {
        label: 'Delete this step',
        color: 'secondary',
        icon: Clear,
        onClick: () => {
          removeStep({
            stepId: step.id,
            doRedirect: {
              currentStepId,
              currentProjectId,
            },
          });
        },
      },
      {
        label: 'Add a step after this step',
        color: 'primary',
        icon: Add,
        onClick: () => {
          updateAddStepDialog({
            $set: {
              open: true,
              showSteps: true,
              showRecipes: false,
              projectId: currentProjectId,
              atStepRank: step.rank + 1,
            },
          });
        },
      },
    ],
  })),
  withStateUpdater('isHovering', false),
  withStyles(styles),
);

export default enhance(StepListItem);
