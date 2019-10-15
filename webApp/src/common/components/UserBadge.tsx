import { Avatar, Button, Dialog, DialogActions, DialogContent, TextField, Typography as T } from '@material-ui/core';
import gql from 'graphql-tag';
import React from 'react';
import { graphql as apolloQuery } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, renderNothing, withProps } from 'recompose';
import Link from 'redux-first-router-link';

import { RFR_LOGOUT } from '../constants';
import { withStateUpdater } from '../utils';

const RetrieveSandboxDialog = compose<any, any>(withStateUpdater('sandboxId', ''))(
  ({ sandboxId, setSandboxId, open, onClose, ...props }) => (
    <Dialog open={open} onClose={onClose}>
      <DialogContent>
        <TextField
          fullWidth
          label="Sandbox ID"
          InputLabelProps={{
            shrink: true,
          }}
          value={sandboxId}
          placeholder="Enter the UUID for the sandbox"
          onChange={(e) => {
            setSandboxId(e.target.value);
          }}
        />
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => {
            if (typeof window === 'object') {
              window.location.href = `/login/debug?role=logged_in_user&aud=postgraphile&user_id=${sandboxId}`;
            }
          }}
        >
          Retrieve
        </Button>
      </DialogActions>
    </Dialog>
  ),
);

const UserBadge: React.FunctionComponent<any> = ({ retrieveSandboxDialog, updateRetrieveSandboxDialog, whoami }) => (
  <div>
    {whoami == null ? (
      <Button component={(props) => <Link to="login" {...props} />}>Login</Button>
    ) : whoami.email === 'sandbox' ? (
      <div style={{ margin: 16 }}>
        <T variant="h6" gutterBottom>
          You are anonymous.
        </T>
        <T gutterBottom>Projects that are not active will be deleted periodically.</T>
        <T variant="caption">
          Your sandbox Id (use it to retrieve the sandbox when the browser cookies are lost or when switching to another
          device):
        </T>
        <T variant="caption">{whoami.id}</T>
        <div style={{ display: 'flex', justifyContent: 'flex-end', margin: '0 -4px 0 -4px' }}>
          <Button
            onClick={() => {
              updateRetrieveSandboxDialog({ open: { $set: true } });
            }}
            style={{ margin: '0 4px' }}
            size="small"
          >
            Retrieve Sandbox
          </Button>
          <Button
            onClick={() => {
              if (typeof window === 'object') {
                window.location.href = '/logout';
              }
            }}
            style={{ margin: '0 4px' }}
            size="small"
          >
            Reset
          </Button>
        </div>
        {retrieveSandboxDialog && (
          <RetrieveSandboxDialog
            open={retrieveSandboxDialog.open}
            onClose={() => {
              updateRetrieveSandboxDialog({ open: { $set: false } });
            }}
          />
        )}
      </div>
    ) : (
      <div
        style={{
          padding: 16,
        }}
      >
        <Avatar
          style={{
            width: '100px',
            height: '100px',
            marginTop: 6,
            marginBottom: 6,
            backgroundSize: 'cover',
            backgroundPosition: 'center',
            backgroundImage: `url(${whoami.profile.picture})`,
            cursor: 'pointer',
          }}
          onClick={() => {
            if (window) {
              window.location.href = '/login';
            }
          }}
        />
        <div style={{ color: '#fff', marginTop: 16 }}>
          <T variant="h4">{whoami.profile.displayName}</T>
          <T variant="body1">{whoami.profile.institution}</T>
        </div>
      </div>
    )}
  </div>
);

const enhance: any = compose(
  // queryRenderer(graphql`
  //   query UserBadgeQuery {
  //     whoami {
  //       id
  //       email
  //       profile {
  //         displayName
  //         picture
  //         institution
  //       }
  //     }
  //   }
  // `),
  connect((state: IReduxState) => ({
    // currentProjectRank: state.location.payload.projectRank,
    currentProjectId: state.app.currentProjectId,
    currentStepId: state.app.currentStepId,
  })),
  apolloQuery(gql`
    query UserBadge {
      whoami {
        id
        email
        profile {
          displayName
          picture
          institution
        }
      }
    }
  `),
  branch((props: any) => props.data.loading, renderNothing),
  branch((props: any) => props.data.whoami == null, renderNothing),
  withStateUpdater('retrieveSandboxDialog', { open: false }),
  withProps(({ data: { whoami } }) => ({
    whoami,
  })),
  // branch(props => props.whoami == null, renderNothing),
);

export default enhance(UserBadge);
