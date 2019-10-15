import { Paper } from '@material-ui/core';
import React from 'react';

import DebugPanel from './DebugPanel';

// TODO: Change this to a real login
const LoginLayout = ({ users }) => (
  <Paper
    elevation={24}
    style={{
      position: 'relative',
      top: '50vh',
      left: '50vw',
      width: 500,
      maxHeight: '90vh',
      overflowY: 'auto',
      transform: 'translate(-50%, -50%)',
      display: 'inline-block',
    }}
  >
    <DebugPanel users={users} />
  </Paper>
);

export default LoginLayout;
