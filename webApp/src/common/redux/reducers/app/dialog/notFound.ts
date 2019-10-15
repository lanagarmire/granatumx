import { CLOSE_NOT_FOUND_DIALOG, OPEN_NOT_FOUND_DIALOG } from '../../../../constants';

export default (state = false, action) => {
  switch (action.type) {
    case OPEN_NOT_FOUND_DIALOG:
      return true;
    case CLOSE_NOT_FOUND_DIALOG:
      return false;
    default:
      return state;
  }
};
