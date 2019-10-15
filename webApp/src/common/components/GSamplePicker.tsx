import { withStyles, TextField, Typography as T } from '@material-ui/core';
import React from 'react';
import { branch, lifecycle, renderComponent, withProps } from 'recompose';
import { compose } from 'redux';
import { queryBackend, withStateUpdater } from '../utils';

const styles: any = {
  // item: {
  //   display: 'block',
  //   verticalAlign: 'baseline !important',
  // },
  // info: {
  //   margin: {
  //     left: '0.2em',
  //   },
  //   fontSize: '0.7em',
  //   opacity: 0.5,
  // },
};

const GSamplePicker = ({
  isLoading,
  data,
  coordinateSelections,
  sampleCoords,
  coords,
}: {
  isLoading: boolean;
  data: Array<{ x: number; y: number; color?: string; [xx: string]: any }>;
  coordinateSelections: Array<{ name: string; id: string }>;
  [xx: string]: any;
}) =>
  isLoading ? (
    <div>Loading ... (GSamplePicker)</div>
  ) : (
    <div>
      <T>Loaded (GSamplePicker)</T>
      <T component={'div' as any}>
        <ul>
          {sampleCoords.map((x) => (
            <li key={JSON.stringify(x)}>{JSON.stringify(x)}</li>
          ))}
        </ul>
      </T>
    </div>
  );

const enhance: any = compose(
  branch(
    ({ gAvailImps }) => gAvailImps == null,
    renderComponent<any>(({ label }) => (
      <TextField label={label} disabled fullWidth value="There are no coordinates. (No available imports)" />
    )),
  ),
  withProps(({ sampleCoords, gAvailImps }) => ({
    sampleCoords: gAvailImps.filter((x) => x.kind === 'sampleCoords'),
  })),
  branch(
    ({ sampleCoords }) => sampleCoords.length === 0,
    renderComponent<any>(({ label, kind }) => (
      <TextField label={label} disabled fullWidth value="There are no coordinates." />
    )),
  ),
  withStateUpdater('isLoading', true),
  withStateUpdater('coords', null),
  lifecycle<{ setIsLoading; setCoords; sampleCoords }, any>({
    async componentDidMount() {
      const coords = (await queryBackend(fetch, {
        audience: '__granatum',
        endpoint: 'getMultipleExports',
        exportIds: [this.props.sampleCoords[0].id],
      }))[0];

      this.props.setCoords(coords);

      this.props.setIsLoading(false);
    },
  }),
  withStyles(styles),
);

export default enhance(GSamplePicker);
