import {
  Button,
  Checkbox,
  FormControlLabel,
  IconButton,
  InputAdornment,
  LinearProgress,
  MenuItem,
  TextField,
  Typography as T,
  withStyles,
} from '@material-ui/core';
import { AttachFile } from '@material-ui/icons';
import filesize from 'filesize';
import React from 'react';
import Dropzone from 'react-dropzone';
import { compose, lifecycle } from 'recompose';

import GSelector from './GSelector';
import GSamplePicker from './GSamplePicker';

const styles = {
  randomizeButton: {
    marginLeft: 16,
  },
};

const ArgInput = ({ section, comp, setState, getState, gAvailImps, classes }) =>
  section === 'imports' ? (
    <GSelector
      kind={comp.kind}
      helperText={comp.description}
      label={comp.label}
      gAvailImps={gAvailImps}
      currentSelection={getState()}
      setCurrentSelection={setState}
    />
  ) : section === 'uploadedFiles' ? (
    <Dropzone
      style={{}}
      onDrop={(accepted) => {
        setState(accepted[0], ['file']);
      }}
      multiple={false}
    >
      <TextField
        fullWidth
        label={`${comp.optional ? '(Optional) ' : ''}${comp.label}`}
        inputProps={{ disabled: true, style: { cursor: 'pointer' } }}
        InputLabelProps={{ shrink: true }}
        InputProps={{
          endAdornment: (
            <InputAdornment position="end">
              <IconButton>
                <AttachFile />
              </IconButton>
            </InputAdornment>
          ),
        }}
        FormHelperTextProps={{
          component: 'div' as any,
          style: { minHeight: 0, marginTop: 0 },
        }}
        helperText={comp.description}
        value={getState(['file']) ? `${getState(['file']).name} (${filesize(getState(['file']).size)})` : ''}
      />
      {getState(['total']) && (
        <LinearProgress
          color={getState(['done']) ? 'primary' : 'secondary'}
          variant="determinate"
          value={(getState(['loaded']) * 100) / getState(['total'])}
        />
      )}
    </Dropzone>
  ) : section === 'args' ? (
    comp.type === 'select' ? (
      <TextField
        onChange={(e) => {
          setState(e.target.value);
        }}
        select
        fullWidth
        InputLabelProps={{
          shrink: true,
        }}
        label={comp.label}
        value={getState() || ''}
        helperText={getState() ? comp.choices.find((x) => x.value === getState()).description : null}
      >
        {comp.choices.map((c) => (
          <MenuItem key={c.value} title={c.description} value={c.value}>
            {c.label || c.value}
          </MenuItem>
        ))}
      </TextField>
    ) : comp.type === 'text' ? (
      <TextField
        fullWidth
        label={`${comp.optional ? '(Optional) ' : ''}${comp.label}`}
        type="text"
        inputProps={{
          maxLength: comp.maxLength,
          pattern: comp.pattern,
        }}
        // eslint-disable-next-line react/jsx-no-duplicate-props
        InputProps={{
          endAdornment: comp.endAdornmentText ? (
            <InputAdornment position="end">{comp.endAdornmentText}</InputAdornment>
          ) : null,
        }}
        value={getState() || ''}
        onChange={(e) => {
          setState(e.target.value);
        }}
        helperText={comp.description}
      />
    ) : comp.type === 'number' ? (
      <TextField
        fullWidth
        label={comp.label}
        type="number"
        inputProps={{
          min: comp.min,
          max: comp.max,
          step: comp.step,
        }}
        // eslint-disable-next-line react/jsx-no-duplicate-props
        InputProps={{
          endAdornment: comp.endAdornmentText ? (
            <InputAdornment position="end">{comp.endAdornmentText}</InputAdornment>
          ) : null,
        }}
        // the state itself is number
        value={getState() != null ? getState().toString() : ''}
        onChange={(e) => {
          setState(+e.target.value);
        }}
        helperText={comp.description}
      />
    ) : comp.type === 'checkbox' ? (
      <div>
        <FormControlLabel
          control={
            <Checkbox
              // checked={console.log('getState() =', getState()) || getState()}
              checked={getState() || false}
              onChange={(e) => {
                console.log('Checkbox:onChange!');
                setState(e.target.checked);
              }}
            />
          }
          label={comp.label}
        />
        {comp.description && <T variant="caption">{comp.description}</T>}
      </div>
    ) : comp.type === 'seed' ? (
      <div style={{ display: 'flex' }}>
        <TextField
          fullWidth
          label={comp.label}
          type="number"
          style={{ flex: 1 }}
          inputProps={{
            min: 0,
            max: 65535,
          }}
          // the state itself is number
          value={(getState() && getState().toString()) || ''}
          onChange={(e) => {
            setState(+e.target.value);
          }}
          helperText={comp.description}
        />
        <Button
          className={classes.randomizeButton}
          onClick={() => {
            setState(Math.floor(Math.random() * Math.floor(65536)));
          }}
        >
          Randomize
        </Button>
      </div>
    ) : comp.type === 'samplePicker' ? (
      <GSamplePicker gAvailImps={gAvailImps} getState={getState()} setState={setState} />
    ) : (
      <T>
        We don&apos;t know how to display an arg with <strong>type: {comp.type}</strong>
      </T>
    )
  ) : (
    <T>
      Unknown <strong>section: {section}</strong>
    </T>
  );

const enhance: any = compose(
  lifecycle<any, any>({
    componentDidMount() {
      const { defaultState, getState, setState } = this.props;
      if (typeof getState() === 'undefined' && typeof defaultState !== 'undefined') {
        setState(defaultState);
      }
    },
  }),
  withStyles(styles),
);

export default enhance(ArgInput);
