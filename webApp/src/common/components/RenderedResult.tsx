import { Table, TableBody, TableCell, TableHead, TableRow, Typography as T, withStyles } from '@material-ui/core';
import React from 'react';
import JSONTree from 'react-json-tree';
import ReactMarkdown from 'react-markdown';
import { CartesianGrid, Scatter, ScatterChart, Tooltip as RTooltip, XAxis, YAxis } from 'recharts';
import { branch, compose, renderNothing } from 'recompose';

// import Heatmap from './Heatmap';
import DataTable from './DataTable';
import ZoomableImage from './ZoomableImage';
import ZScatter from './ZScatter';

const styles: any = (theme) => ({
  source: {
    margin: [theme.spacing.unit * 2, 0],
    fontFamily: 'monospace',
    fontSize: '1em',
    whiteSpace: 'pre-wrap',
    wordBreak: 'break-all',
    wordWrap: 'break-word',
  },
  png: {
    textAlign: 'center',
  },
  sourceEllipsis: {
    opacity: 0.3,
  },
  sourceOmission: {
    margin: [10, 30],
    opacity: 0.3,
  },
  divCenter: {
    marginLeft: 'auto',
    marginRight: 'auto',
  },
});

const DEFAULT_HEIGHT = 650;
const DEFAULT_WIDTH = 750;

const RenderedResult: React.FunctionComponent<any> = ({ type, data, classes, ...props }) =>
  type === 'iframe' ? (
    <div style={{ textAlign: 'center' }}>
      <iframe
        style={{ display: 'block', margin: 'auto' }}
        title={props.title}
        frameBorder="0"
        width={props.width || DEFAULT_WIDTH}
        height={props.height || DEFAULT_HEIGHT}
        srcDoc={data}
      />
      <T variant="caption">{props.description}</T>
    </div>
  ) : type === 'scatterplot' ? (
    <div style={{ textAlign: 'center' }}>
      <ScatterChart
        width={props.width || DEFAULT_WIDTH}
        height={props.height || DEFAULT_HEIGHT}
        margin={{ top: 20, right: 20, bottom: 10, left: 10 }}
      >
        <CartesianGrid strokeDasharray="3 3" />
        <XAxis dataKey="x" name={data.dimNames[0]} type="number" />
        <YAxis dataKey="y" name={data.dimNames[1]} type="number" />
        <Scatter
          isAnimationActive={false}
          data={Object.entries(data.coords).map(([k, v]) => ({ id: k, x: v[0], y: v[1] }))}
          fill="#8884d8"
        />
        <RTooltip isAnimationActive={false} cursor={{ strokeDasharray: '3 3' }} />
      </ScatterChart>
      <T variant="caption">{props.description}</T>
    </div>
  ) : type === 'png' ? (
    <ZoomableImage
      type="png"
      width={props.width || DEFAULT_WIDTH}
      height={props.height || DEFAULT_HEIGHT}
      data={data}
      description={props.description}
      className={classes.png}
    />
  ) : type === 'code' ? (
    <div className={classes.source}>{data}</div>
  ) : type === 'codeExcerpt' ? (
    <div className={classes.source}>
      {data.head}
      <span className={classes.sourceEllipsis}>...</span>
      <div className={classes.sourceOmission}>(omitting {data.numBytesOmitted} bytes)</div>
      <span className={classes.sourceEllipsis}>...</span>
      {data.tail}
    </div>
  ) : type === 'table' ? (
    <DataTable {...{ data, ...props }} />
  ) : type === 'compactTable' ? (
    <Table className={classes.table}>
      <TableHead>
        <TableRow>
          {data.cols.map((c) => (
            <TableCell key={c}>{c}</TableCell>
          ))}
        </TableRow>
      </TableHead>
      <TableBody>
        {data.data.map((r, i) => (
          <TableRow key={r._id || i}>
            {data.cols.map((c) => (
              <TableCell key={c}>{r[c]}</TableCell>
            ))}
          </TableRow>
        ))}
      </TableBody>
    </Table>
  ) : type === 'heatmap' ? (
    <div>
      {/* <Heatmap data={data} width={props.width || DEFAULT_WIDTH} height={props.height || DEFAULT_HEIGHT} /> */}
    </div>
  ) : type === 'markdown' ? (
    <T component={'div' as any}>
      <ReactMarkdown source={data} />
    </T>
  ) : type === 'text' ? (
    <T>
      {typeof props.label === 'string' ? `${props.label}: ` : ''}
      {typeof data === 'string' ? data : `[ERROR: typeof data === ${typeof data}]`}
    </T>
  ) : (
    <JSONTree data={data} />
  );

const enhance = compose(
  branch(({ data }) => data == null, renderNothing),
  withStyles(styles),
);

export default enhance(RenderedResult);
