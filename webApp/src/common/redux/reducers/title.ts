import {
  RFR_CREATE_PROJECT,
  RFR_NO_STEPS,
  RFR_NO_STEPS_THUNK_END,
  RFR_PROJECT_STEP,
  RFR_PROJECT_STEP_THUNK_END,
} from '../../constants';

export default (state = 'GranatumX', action) => {
  switch (action.type) {
    case RFR_CREATE_PROJECT:
      return 'GranatumX - Create a project';
    case RFR_NO_STEPS: {
      const { projectId } = action.payload;
      return `GranatumX - [Project: ${projectId}] - No step selected`;
    }
    case RFR_PROJECT_STEP: {
      const { stepId } = action.payload;
      return `GranatumX - [Step ${stepId}]`;
    }
    case RFR_NO_STEPS_THUNK_END: {
      const { project } = action.payload;
      return `GranatumX - [Project: ${project.name}] - No step selected`;
    }
    case RFR_PROJECT_STEP_THUNK_END: {
      const { project, step } = action.payload;
      return `GranatumX - [Project: ${project.name}] - [Step ${step.rank + 1}. ${step.gboxByGbox.meta.title}]`;
    }
    default:
      return state;
  }
};
