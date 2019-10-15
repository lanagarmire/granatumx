import { Chip, withWidth } from '@material-ui/core';
import { withStyles } from '@material-ui/core/es';
import { KeyboardArrowRight } from '@material-ui/icons';
import React from 'react';
import { compose } from 'recompose';

const Pipeline: React.FunctionComponent<any> = ({ steps, ...props }) => (
  <div style={{ margin: '8px -8px', display: 'flex', alignItems: 'center', flexFlow: 'row wrap' }} {...props}>
    {steps.map((step, i) => (
      <span key={step.id}>
        <Chip style={{ margin: 8 }} label={step.title} />
        {i < steps.length - 1 && <KeyboardArrowRight style={{ margin: -8 }} />}
      </span>
    ))}
  </div>
);

const enhance: any = compose();

export default enhance(Pipeline);
