import update from 'immutability-helper';

import { UPDATE_ADD_STEP_DIALOG } from '../../../../constants';

export default (
  state = { open: false, showSteps: true, showRecipes: true, projectId: undefined, atStepRank: undefined },
  action,
) => {
  switch (action.type) {
    case UPDATE_ADD_STEP_DIALOG:
      return update(state, action.payload);
    default:
      return state;
  }
};
