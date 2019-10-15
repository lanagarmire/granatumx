import { Typography as T } from '@material-ui/core';
import { NoteAdd } from '@material-ui/icons';
import React from 'react';
import { connect } from 'react-redux';
import { NOT_FOUND } from 'redux-first-router';

import { RFR_CREATE_PROJECT, RFR_NO_STEPS, RFR_PROJECT_REPORT, RFR_PROJECT_ROOT, RFR_PROJECT_STEP } from '../constants';
import CreateProjectPage from './CreateProjectPage';
import NotFoundDialog from './NotFoundDialog';
import ProjectStepLayout from './ProjectStepLayout';
import ReportLayout from './ReportLayout';

const LayoutChooser = ({ type }) =>
  type === RFR_NO_STEPS ? (
    <ProjectStepLayout />
  ) : type === RFR_PROJECT_ROOT ? (
    <div style={{ marginTop: 200 }}>
      <NoteAdd style={{ fontSize: '64px' }} />
      <T variant="h2" align="center">
        No step is selected. Please select or create steps from the left panel.
      </T>
    </div>
  ) : type === RFR_CREATE_PROJECT ? (
    <div
      style={{
        position: 'relative',
        top: '50vh',
        left: '50vw',
        width: '450px',
        transform: 'translate(-50%, -50%)',
        display: 'inline-block',
      }}
    >
      <CreateProjectPage />
    </div>
  ) : type === RFR_PROJECT_REPORT ? (
    <ReportLayout />
  ) : type === RFR_PROJECT_STEP ? (
    <ProjectStepLayout />
  ) : type === NOT_FOUND ? (
    <NotFoundDialog />
  ) : null;

export default connect((state: IReduxState) => ({ type: state.location.type }))(LayoutChooser);
