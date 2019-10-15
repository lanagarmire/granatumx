import { MenuItem, TextField, withStyles } from '@material-ui/core';
import React from 'react';
import { branch, lifecycle, renderComponent, withProps } from 'recompose';
import { compose } from 'redux';

const styles: any = {
  item: {
    display: 'block',
    verticalAlign: 'baseline !important',
  },
  info: {
    margin: {
      left: '0.2em',
    },
    fontSize: '0.7em',
    opacity: 0.5,
  },
};

const GSelector = ({ helperText, label, imports, currentSelection, setCurrentSelection, classes }) => (
  <TextField
    label={label}
    fullWidth
    helperText={helperText}
    select
    value={currentSelection || ''}
    onChange={(e) => setCurrentSelection(e.target.value)}
  >
    {imports.map((x) => (
      <MenuItem key={x.id} value={x.id} classes={{ root: classes.item }}>
        {x.extractFrom}
        <span className={classes.info}>
          (from Step {x.stepRank + 1}: {x.stepTitle})
        </span>
      </MenuItem>
    ))}
  </TextField>
);

const enhance: any = compose(
  branch(
    ({ gAvailImps }) => gAvailImps == null,
    renderComponent<any>(({ label }) => (
      <TextField label={label} disabled fullWidth value="Available imports not found." />
    )),
  ),
  withProps(({ kind, gAvailImps }) => ({
    imports: gAvailImps.filter((x) => x.kind === kind),
  })),
  branch(
    ({ imports }) => imports.length === 0,
    renderComponent<any>(({ label, kind }) => (
      <TextField label={label} disabled fullWidth value={`There are no available ${kind} to use.`} />
    )),
  ),
  lifecycle<{ imports; currentSelection; setCurrentSelection }, any>({
    componentDidMount() {
      const { imports, currentSelection, setCurrentSelection } = this.props;
      if (currentSelection == null) {
        setCurrentSelection(imports[imports.length - 1].id);
      }
    },
  }),
  withStyles(styles),
);

export default enhance(GSelector);
