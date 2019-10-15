import { List, ListItem, ListItemText } from '@material-ui/core';
import React from 'react';

const DebugPanel = ({ users }) => (
  <div style={{ overflowY: 'auto' }}>
    <List>
      {users.map((user) => (
        <ListItem
          button
          component={(props) => (
            // eslint-disable-next-line jsx-a11y/anchor-has-content
            <a href={`/login/debug?role=logged_in_user&aud=postgraphile&user_id=${user.id}`} {...props as any} />
          )}
          key={user.id}
        >
          <ListItemText primary={user.id} secondary={`${user.email} - ${user.display_name} - ${user.institution}`} />
        </ListItem>
      ))}
      <ListItem
        component={(props) => (
          // eslint-disable-next-line jsx-a11y/anchor-has-content
          <a href="/logout" {...props as any} />
        )}
        button
        key="logout"
      >
        <ListItemText primary="Logout" />
      </ListItem>
    </List>
  </div>
);

export default DebugPanel;
