import { Map } from 'immutable';

import { COMPONENT_SET_STATE } from '../../../constants';

export default (state: Map<string, any> = Map(), action) => {
  switch (action.type) {
    case COMPONENT_SET_STATE:
      return state.setIn(action.payload.path, action.payload.newState);
    default:
      return state;
  }
};
