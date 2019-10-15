import { RFR_NO_STEPS_THUNK_END, RFR_PROJECT_STEP_THUNK_END } from '../../../constants';

export default (state = null, action) => {
  switch (action.type) {
    case RFR_PROJECT_STEP_THUNK_END:
      return action.payload.step.id;
    case RFR_NO_STEPS_THUNK_END:
      return null;
    default:
      return state;
  }
};
