import uuid from 'uuid';
import { DISMISS_GLOBAL_DIALOG, GLOBAL_DIALOG, GLOBAL_ERROR } from '../../constants';

export default (state: IStateGlobalDialog[] = [], action) => {
  switch (action.type) {
    case GLOBAL_ERROR:
      return [...state, { type: 'error', id: uuid.v4(), ...action.payload }];
    case GLOBAL_DIALOG:
      return [...state, { type: null, id: uuid.v4(), ...action.payload }];
    case DISMISS_GLOBAL_DIALOG: {
      const { id } = action.payload;
      return state.filter((e) => id !== e.id);
    }
    default:
      return state;
  }
};
