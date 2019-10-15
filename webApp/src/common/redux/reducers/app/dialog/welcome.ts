import update from 'immutability-helper';

import { UPDATE_WELCOME_DIALOG } from '../../../../constants';

export default (state = { open: false }, action) => {
  switch (action.type) {
    case UPDATE_WELCOME_DIALOG:
      return update(state, action.payload);
    default:
      return state;
  }
};
