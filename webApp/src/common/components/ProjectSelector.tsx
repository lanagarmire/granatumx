import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  TextField,
  Typography as T,
  withStyles,
} from '@material-ui/core';
import {
  Accessibility,
  Assignment,
  ClearAll,
  CloudDownload,
  Delete,
  FiberNew,
  Info,
  Menu as MenuIcon,
} from '@material-ui/icons';
import gql from 'graphql-tag';
import _ from 'lodash';
import React from 'react';
import { graphql } from 'react-apollo';
import ReactMarkdown from 'react-markdown';
import { connect } from 'react-redux';
import { branch, compose, renderComponent, renderNothing, withProps, withState } from 'recompose';

import { DESCRIPTION_LIMIT, NAME_LIMIT } from '../constants';
import { withStateUpdater } from '../utils';
import CreateProjectPage from './CreateProjectPage';
import SavePipelineDialog from './SavePipelineDialog';
import actionCreators from '../redux/actionCreators';

const NoProjects = () => (
  <div style={{ padding: 16 }}>
    <T>No projects were found.</T>
  </div>
);

const styles = {
  toolbar: {
    display: 'flex',
    justifyContent: 'flex-end',
    marginTop: 16,
    minHeight: 36,
  },
  toolbarButton: {
    width: 36,
    height: 36,
    '&:hover': {
      background: '#f7f7f7',
    },
    padding: 0,
    marginLeft: 16,
  },
  createProjectDialogPaper: {
    width: 450,
  },
};

const EditingDialog: React.ComponentClass<any> = compose<any, any>(
  withState('name', 'setName', ({ name }) => name),
  withState('description', 'setDescription', ({ description }) => description),
  withProps(({ name }) => ({
    nameError: !name
      ? 'Cannot be empty!'
      : name.length > NAME_LIMIT
      ? `Too long (limit ${NAME_LIMIT} characters)!`
      : '',
  })),
  withProps(({ description }) => ({
    descriptionError: description.length > DESCRIPTION_LIMIT ? `Too long (limit ${DESCRIPTION_LIMIT} characters)!` : '',
  })),
  withStyles({
    projectDetailsDialog: {
      width: '700px',
      maxWidth: 'none',
    },
  }),
)(({ submit, open, onClose, classes, name, setName, nameError, description, setDescription, descriptionError }) => (
  <Dialog classes={{ paper: classes.projectDetailsDialog }} open={open}>
    <DialogTitle>Edit project info</DialogTitle>
    <DialogContent>
      <TextField
        fullWidth
        error={!!nameError}
        helperText={nameError}
        label="Project name"
        InputLabelProps={{
          shrink: true,
        }}
        value={name}
        placeholder="My project"
        onChange={(e) => {
          setName(e.target.value);
        }}
      />
    </DialogContent>
    <DialogContent>
      <TextField
        fullWidth
        multiline
        error={!!descriptionError}
        helperText={descriptionError}
        rows="30"
        label="Project description"
        InputLabelProps={{
          shrink: true,
        }}
        value={description}
        placeholder={
          'Write a brief description of the project. The text is parsed as a markdown. You can write down\n' +
          '\n' +
          '  * The objectives\n' +
          '  * Experiment design\n' +
          '  * The source of the data\n' +
          '  * The expected results\n' +
          '  * The planned analysis steps\n' +
          '\n' +
          'You can add more comments.'
        }
        onChange={(e) => {
          setDescription(e.target.value);
        }}
        InputProps={{
          style: {
            fontFamily: 'monospace',
          },
        }}
      />
    </DialogContent>
    <DialogActions>
      <Button
        onClick={() => {
          onClose();
        }}
      >
        Cancel
      </Button>
      <Button
        color="primary"
        onClick={() => {
          submit({ name, description });
        }}
      >
        Save
      </Button>
    </DialogActions>
  </Dialog>
));

const ProjectSelector: React.FunctionComponent<any> = ({
  classes,
  createProjectDialog,
  currentProject,
  currentProjectId,
  detailsDialog,
  editingDialog,
  menuAnchorEl,
  menuItems,
  modifyProject,
  projects,
  savePipelineDialog,
  setCreateProjectDialog,
  setDetailsDialog,
  setEditingDialog,
  setMenuAnchorEl,
  rfrProjectRoot,
  updateSavePipelineDialog,
}) => (
  <div style={{ padding: 16, paddingTop: 0 }}>
    <div style={{ display: 'flex' }}>
      <TextField
        style={{ flex: 1 }}
        onChange={(e) => {
          rfrProjectRoot({ projectId: e.target.value });
        }}
        id="project"
        select
        fullWidth
        InputLabelProps={{
          shrink: true,
        }}
        label="Project"
        value={currentProjectId}
      >
        {projects.map((project) => (
          <MenuItem key={project.id} value={project.id}>
            {project.name}
          </MenuItem>
        ))}
      </TextField>
      <IconButton
        style={{ alignSelf: 'flex-end' }}
        className={classes.toolbarButton}
        onClick={(e) => {
          setMenuAnchorEl(e.target);
        }}
      >
        <MenuIcon />
      </IconButton>
    </div>
    <Menu
      anchorEl={menuAnchorEl}
      open={Boolean(menuAnchorEl)}
      onClose={() => {
        setMenuAnchorEl(null);
      }}
    >
      {menuItems.map((m) => (
        <MenuItem
          key={m.title}
          onClick={(e) => {
            setMenuAnchorEl(null);
            m.onClick(e);
          }}
        >
          <ListItemIcon>{m.icon}</ListItemIcon>
          <ListItemText primary={m.title} />
        </MenuItem>
      ))}
    </Menu>
    {savePipelineDialog.open && (
      <SavePipelineDialog
        open={savePipelineDialog.open}
        onClose={() => {
          updateSavePipelineDialog({ open: { $set: false } });
        }}
      />
    )}
    {editingDialog && (
      <EditingDialog
        open={editingDialog}
        onClose={() => {
          setEditingDialog(false);
          setDetailsDialog(true);
        }}
        submit={({ name, description }) => {
          modifyProject({ id: currentProjectId, name, description });
          setEditingDialog(false);
          setDetailsDialog(true);
        }}
        name={currentProject.name}
        description={currentProject.description}
      />
    )}
    <Dialog
      classes={{ paper: classes.projectDetailsDialog }}
      open={detailsDialog}
      onClose={() => setDetailsDialog(false)}
    >
      <DialogTitle>{currentProject.name}</DialogTitle>
      <DialogContent>
        <T component={'div' as any}>
          <ReactMarkdown source={currentProject.description} />
        </T>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => {
            setDetailsDialog(false);
            setEditingDialog(true);
          }}
        >
          Edit
        </Button>
        <Button onClick={() => setDetailsDialog(false)}>OK</Button>
      </DialogActions>
    </Dialog>
    <Dialog
      open={createProjectDialog}
      classes={{ paper: classes.createProjectDialogPaper }}
      onClose={() => setCreateProjectDialog(false)}
    >
      <CreateProjectPage asDialog onClose={() => setCreateProjectDialog(false)} />
    </Dialog>
  </div>
);

const enhance = compose(
  connect((state: IReduxState) => ({
    currentProjectId: state.app.currentProjectId,
  })),
  branch((props: any) => props.currentProjectId == null, renderNothing),
  graphql(
    gql`
      query ProjectSelector {
        whoami {
          id
          projectsByOwnerId {
            nodes {
              id
              rank
              name
              description
            }
          }
        }
      }
    `,
  ),
  branch((props: any) => props.data.loading, renderNothing),
  branch((props: any) => props.data.whoami == null, renderNothing),
  withProps(({ data }) => ({
    projects: _(data.whoami.projectsByOwnerId.nodes)
      .orderBy('rank')
      .value(),
  })),
  withProps(({ projects, currentProjectId }) => ({
    currentProject: projects.find((x) => x.id === currentProjectId),
  })),
  branch(({ currentProject }) => currentProject == null, renderNothing),
  branch(({ projects }) => projects.length === 0, renderComponent(NoProjects)),
  connect(
    null,
    {
      rfrProjectRoot: actionCreators.rfrProjectRoot,
      rfrNoSteps: actionCreators.rfrNoSteps,
      deleteProject: actionCreators.deleteProject,
      modifyProject: actionCreators.modifyProject,
      removeAllSteps: actionCreators.removeAllSteps,
      rfrProjectReport: actionCreators.rfrProjectReport,
      updateWelcomeDialog: actionCreators.updateWelcomeDialog,
    },
  ),
  withStateUpdater('detailsDialog', false),
  withStateUpdater('editingDialog', false),
  withStateUpdater('createProjectDialog', false),
  withStateUpdater('savePipelineDialog', { open: false }),
  withStateUpdater('menuAnchorEl', null),
  withProps(
    ({
      rfrProjectReport,
      removeAllSteps,
      rfrNoSteps,
      deleteProject,
      currentProjectId,
      setMenuAnchorEl,
      setCreateProjectDialog,
      updateSavePipelineDialog,
      setDetailsDialog,
      updateWelcomeDialog,
    }) => ({
      menuItems: [
        {
          title: 'Show the welcome message',
          icon: <Accessibility />,
          onClick: () => {
            updateWelcomeDialog({ open: { $set: true } });
          },
        },
        {
          title: 'Project report',
          icon: <Assignment />,
          onClick: () => {
            rfrProjectReport({ projectId: currentProjectId });
          },
        },
        {
          title: 'Delete this project',
          icon: <Delete />,
          onClick: () => {
            deleteProject({ projectId: currentProjectId });
          },
        },
        {
          title: 'Create a new project',
          icon: <FiberNew />,
          onClick: () => {
            setCreateProjectDialog(true);
          },
        },
        {
          title: 'Show project description',
          icon: <Info />,
          onClick: () => {
            setDetailsDialog(true);
          },
        },
        {
          title: 'Save current pipeline as a recipe',
          icon: <CloudDownload />,
          onClick: () => {
            updateSavePipelineDialog({ open: { $set: true } });
          },
        },
        {
          title: 'Remove all steps',
          icon: <ClearAll />,
          onClick: () => {
            setMenuAnchorEl(null);
            if (confirm('Are you sure? All data in this project will be deleted.')) {
              removeAllSteps({ projectId: currentProjectId });
              rfrNoSteps({ projectId: currentProjectId });
            }
          },
        },
      ],
    }),
  ),
  withStyles(styles),
);

export default enhance(ProjectSelector);
