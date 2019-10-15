import { Button, DialogActions, DialogContent, DialogTitle, TextField, withStyles } from '@material-ui/core';
import React from 'react';
import { connect } from 'react-redux';
import { compose, defaultProps, withHandlers, withProps, withState } from 'recompose';

import { CREATE_PROJECT, DESCRIPTION_LIMIT, NAME_LIMIT } from '../constants';
import actionCreators from '../redux/actionCreators';

const styles = {
  buttonRow: {
    display: 'flex',
    margin: '4px -4px',
  },
  button: {
    margin: 4,
  },
};

const CreateProjectPage: React.FunctionComponent<any> = ({
  classes,
  submit,
  onClose,
  name,
  setName,
  nameError,
  description,
  setDescription,
  descriptionError,
  asDialog,
}) => (
  <>
    <DialogTitle>Create a new project</DialogTitle>
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
        rows="10"
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
      {asDialog && <Button onClick={onClose}>Cancel</Button>}
      <Button
        variant={asDialog ? 'text' : 'contained'}
        color="primary"
        onClick={() => {
          onClose();
          submit();
        }}
      >
        Create
      </Button>
    </DialogActions>
  </>
);

const enhance: any = compose(
  withState('name', 'setName', 'My project'),
  withState('description', 'setDescription', ''),
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
  connect(
    null,
    {
      createProject: actionCreators.createProject,
    },
  ),
  withHandlers({
    submit: ({ nameError, name, descriptionError, description, createProject }) => () => {
      if (nameError || descriptionError) {
        return;
      }

      createProject({ name, description });
    },
  }),
  defaultProps({
    onClose: () => null,
    asDialog: false,
  }),
  withStyles(styles),
);

export default enhance(CreateProjectPage);
