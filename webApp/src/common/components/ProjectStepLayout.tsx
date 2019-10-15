/**
 * TODO: IDEA: In the upload stage, let the user tag some popular tags (like "id", "batch") to the metadata
 * TODO: IDEA: Although steps are always linearly ordered, in the end there'd be a dependency tree. The users can pick particular branches
 * TODO: The SampleMetaSelector/GeneMetaSelector should also have the ability to ignore the recommended tags
 * TODO: If the meta column is bool, give them an option to select inversed (for e.g. outliers)
 * The default double-panel layout seen in most parts of the app.
 *
 * The reason this is isolated from the <App /> component is because some
 * pages might have an alternative layout, like "login", "about", or "welcome"
 *
 */

import {
  Button,
  Divider,
  Drawer,
  ExpansionPanel,
  ExpansionPanelDetails,
  ExpansionPanelSummary,
  Tooltip,
  Typography as T,
  withStyles,
  withWidth,
} from '@material-ui/core';
import cn from 'classnames';
import React from 'react';
import { connect } from 'react-redux';
import { compose } from 'redux';

import {
  ChevronLeft,
  ChevronRight,
  ExpandMore as ExpandMoreIcon,
  Menu as MenuIcon,
  ViewList,
} from '@material-ui/icons';
import { withProps } from 'recompose';
import actionCreators from '../redux/actionCreators';
import { withStateUpdater } from '../utils';
import AddStepDialog from './AddStepDialog';
import DataManager from './DataManager';
import GlobalDialogs from './GlobalDialogs';
import ProjectSelector from './ProjectSelector';
import StepDisplayer from './StepDisplayer';
import StepJustFinished from './StepJustFinished';
import StepList from './StepList';
import UserBadge from './UserBadge';
import WelcomeDialog from './WelcomeDialog';
import Logo from './Logo';

const LEFT_DRAWER_WIDTH = 300;
const RIGHT_DRAWER_WIDTH = 300;
const MIDDLE_AREA_WIDTH = 800;

const styles: any = {
  drawer: {
    backgroundSize: 'initial',
    backgroundPosition: 'top',
    display: 'flex',
  },
  drawerLeft: {
    width: LEFT_DRAWER_WIDTH,
  },
  drawerRight: {
    width: RIGHT_DRAWER_WIDTH,
  },
  main: {
    height: '100vh',
    overflow: 'auto',
    backgroundColor: '#fafafa',
    position: 'relative',
  },
  mainRetractLeft: {
    marginLeft: LEFT_DRAWER_WIDTH,
  },
  mainRetractRight: {
    marginRight: RIGHT_DRAWER_WIDTH,
  },
  mobileButton: {
    flex: 1,
  },
  mobileExpansion: {
    boxShadow: 'none',
  },
  mobileUserBadge: {
    padding: 0,
  },
};

const ProjectStepLayout = ({
  setLeftPanelOpen,
  leftPanelOpen,
  setRightPanelOpen,
  rightPanelOpen,
  currentProjectId,
  currentStepId,
  mobileActions,
  isMobile,
  classes,
}) => (
  <div>
    <Drawer
      anchor="left"
      variant={isMobile ? 'temporary' : 'permanent'}
      classes={{
        paper: [classes.drawer, classes.drawerLeft].join(' '),
      }}
      open={!isMobile || leftPanelOpen}
      onClose={() => {
        setLeftPanelOpen(false);
      }}
    >
      {isMobile ? (
        <ExpansionPanel className={classes.mobileExpansion}>
          <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />}>
            <Logo />
          </ExpansionPanelSummary>
          <ExpansionPanelDetails className={classes.mobileUserBadge}>
            <UserBadge />
          </ExpansionPanelDetails>
        </ExpansionPanel>
      ) : (
        <>
          <Logo />
          <UserBadge />
        </>
      )}
      <ProjectSelector />
      <Divider />
      <StepList />
    </Drawer>
    <main className={cn(classes.main, isMobile || classes.mainRetractLeft, isMobile || classes.mainRetractRight)}>
      {isMobile && (
        <div style={{ display: 'flex', height: 64 }}>
          {mobileActions.map((m) => (
            <Tooltip key={m.label} title={m.label}>
              <Button onClick={m.onClick} className={classes.mobileButton}>
                {React.createElement(m.icon)}
              </Button>
            </Tooltip>
          ))}
        </div>
      )}
      {currentProjectId == null || currentStepId == null ? (
        <div
          style={{
            position: 'absolute',
            top: '50%',
            left: '50%',
            transform: 'translate(-50%, -50%)',
          }}
        >
          <div style={{ margin: '32px', textAlign: 'center' }}>
            <img
              style={{
                width: 300,
                marginBottom: 32,
                userSelect: 'none',
                filter: 'grayscale(100%) contrast(30%) brightness(150%)',
              }}
              alt="granatumx_logo_zzz.svg"
              src="/granatumx_logo_zzz.svg"
            />
            <T variant="caption">No step is selected. Select or create a step from the left panel.</T>
          </div>
        </div>
      ) : (
        <div
          style={{
            width: isMobile ? 'calc(100% - 32px)' : MIDDLE_AREA_WIDTH,
            margin: '0 auto',
          }}
        >
          <StepDisplayer {...{ currentProjectId, currentStepId }} />
        </div>
      )}
    </main>
    <Drawer
      classes={{
        paper: [classes.drawer, classes.drawerRight].join(' '),
      }}
      anchor="right"
      variant={isMobile ? 'temporary' : 'permanent'}
      open={!isMobile || rightPanelOpen}
      onClose={() => {
        setRightPanelOpen(false);
      }}
    >
      <div style={{ flex: 1 }}>
        <DataManager />
      </div>
    </Drawer>
    <GlobalDialogs />
    <AddStepDialog />
    <StepJustFinished />
    <WelcomeDialog />
  </div>
);

const enhance = compose(
  connect(
    (state: IReduxState) => ({
      currentStepId: state.app.currentStepId,
      currentProjectId: state.app.currentProjectId,
    }),
    {
      nextStep: actionCreators.nextStep,
      prevStep: actionCreators.prevStep,
    },
  ),
  withStateUpdater('leftPanelOpen', false),
  withStateUpdater('rightPanelOpen', false),
  withProps(({ setLeftPanelOpen, setRightPanelOpen, prevStep, nextStep }) => ({
    mobileActions: [
      {
        label: 'Open step panel',
        icon: MenuIcon,
        onClick: () => {
          setLeftPanelOpen(true);
        },
      },
      {
        label: 'Previous step',
        icon: ChevronLeft,
        onClick: () => {
          prevStep({ reportAtTheEnd: true });
        },
      },
      {
        label: 'Next step',
        icon: ChevronRight,
        onClick: () => {
          nextStep({ reportAtTheEnd: true });
        },
      },
      {
        label: 'Open project data panel',
        icon: ViewList,
        onClick: () => {
          setRightPanelOpen(true);
        },
      },
    ],
  })),
  withWidth(),
  withProps(({ width }) => ({
    isMobile: width === 'xs' || width === 'sm' || width === 'md',
  })),
  withStyles(styles),
);

export default enhance(ProjectStepLayout);
