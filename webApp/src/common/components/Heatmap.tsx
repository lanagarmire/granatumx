import Canvas from 'isomorphic-canvas';
import React from 'react';
import { compose, withState } from 'recompose';

const Heatmap: React.FunctionComponent<any> = ({ data, width, height, zoom, setZoom, offsets, setOffsets }) => {
  const ncols = data.sampleIds.length;
  const nrows = data.geneIds.length;

  const canvas = Canvas(ncols, nrows);
  const ctx = canvas.getContext('2d');
  const imgData = ctx.createImageData(ncols, nrows);

  for (let i = 0; i < nrows; i++) {
    for (let j = 0; j < ncols; j++) {
      imgData.data[i * ncols * 4 + j * 4] = Math.floor(10 * data.matrix[i][j]);
      imgData.data[i * ncols * 4 + j * 4 + 1] = Math.floor(10 * data.matrix[i][j]);
      imgData.data[i * ncols * 4 + j * 4 + 2] = Math.floor(10 * data.matrix[i][j]);
      imgData.data[i * ncols * 4 + j * 4 + 3] = 255;
    }
  }

  ctx.putImageData(imgData, 0, 0);
  const href = canvas.toDataURL('image/png');

  // const zoomFactor = (0.5 * width) / ncols;

  return (
    <div style={{ position: 'relative', textAlign: 'center', overflow: 'auto' }}>
      <svg width={width} height={height} viewBox="0 0 100 100" preserveAspectRatio="none">
        <image
          imageRendering="pixelated"
          width={ncols}
          height={nrows}
          style={{ transform: `scale(${1})` }}
          xlinkHref={href}
          onMouseMove={(e) => {
            setOffsets({ x: e.nativeEvent.offsetX, y: e.nativeEvent.offsetY });
          }}
        />
      </svg>
    </div>
  );
};

const enhance = compose(
  withState('zoom', 'setZoom', null),
  withState('offsets', 'setOffsets', { x: null, y: null }),
);

export default enhance(Heatmap);
