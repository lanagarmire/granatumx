import React from 'react';
import { compose } from 'recompose';

const Logo = () => (
  <div style={{ margin: 'auto', marginTop: 16, marginBottom: 0 }}>
    <img
      style={{
        width: 200,
        userSelect: 'none',
      }}
      alt="granatumx_logo_with_text.svg"
      src="/granatumx_logo_with_text.svg"
    />
  </div>
);

const enhance = compose();

export default enhance(Logo);
