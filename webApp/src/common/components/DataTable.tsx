import { Tooltip, Typography as T, withStyles } from '@material-ui/core';

import { IntegratedSorting, SortingState, DataTypeProvider } from '@devexpress/dx-react-grid';
import {
  Grid,
  TableColumnReordering,
  TableColumnResizing,
  TableHeaderRow,
  VirtualTable,
} from '@devexpress/dx-react-grid-material-ui';

import React from 'react';
import { compose, withHandlers, withProps, withState } from 'recompose';
import _ from 'lodash';
import { range0 } from '../utils';

const styles: any = {
  table: {
    marginBottom: 32,
  },
  tableContainer: {
    overflowX: 'auto',
  },
  caption: {
    color: 'rgba(0, 0, 0, 0.54)',
    textAlign: 'center',
    marginTop: 16,
  },
};

const TooltippedCell: React.FunctionComponent<any> = ({ value }) => (
  <Tooltip title={value != null ? value.toString() : ''}>
    <div>{value}</div>
  </Tooltip>
);

const DataTable: React.FunctionComponent<any> = ({
  numPages,
  prevPage,
  nextPage,
  pageSize,
  classes,
  columns,
  rows,
  title,
  caption,
  page,
  colNames,
  setPage,
  ...props
}) => (
  <div>
    {title != null && <T variant="h6">{title}</T>}
    <div className={classes.tableContainer}>
      {console.log(colNames)}
      <Grid rows={rows} columns={columns}>
        <DataTypeProvider for={colNames} formatterComponent={TooltippedCell} />
        <SortingState />
        <IntegratedSorting />
        <VirtualTable />
        <TableHeaderRow showSortingControls />
      </Grid>
    </div>
    {caption != null && (
      <T variant="caption" className={classes.caption}>
        {caption}
      </T>
    )}
  </div>
);

const enhance: any = compose(
  withProps(({ data: { title, columns, data, pageSize, caption } }) => ({
    title,
    caption,
    columns:
      columns == null
        ? range0(data[0].length).map((x) => ({ name: x.toString(), title: x.toString() }))
        : columns
            .map((x, i) => {
              if (typeof x === 'object') {
                return x;
              } else {
                return { name: x.toString(), title: x.toString() };
              }
            })
            .map((x) => ({ ...x, columnName: x.name, width: 750 / columns.length })),
  })),
  withProps(({ data: { data }, columns }) => ({
    rows: data.map((r) => _.fromPairs(r.map((c, i) => [columns[i].name, c]))),
    colNames: columns.map((x) => x.name),
  })),
  withProps(({ data: { data }, pageSize }) => ({
    numPages: data.length === 0 ? 0 : Math.floor((data.length - 1) / pageSize) + 1,
  })),
  withState('page', 'setPage', 0),
  withHandlers({
    prevPage: ({ page, setPage }) => () => {
      if (page > 0) {
        setPage(page - 1);
      }
    },
    nextPage: ({ numPages, page, setPage }) => () => {
      if (page < numPages - 1) {
        setPage(page + 1);
      }
    },
  }),
  withStyles(styles),
);

export default enhance(DataTable);
