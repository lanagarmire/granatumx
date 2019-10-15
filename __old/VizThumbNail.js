import React from 'react';
import {
  Legend,
  LineChart,
  Line,
  ScatterChart,
  Scatter,
  CartesianGrid,
  XAxis,
  YAxis,
  ResponsiveContainer,
  Brush,
  Tooltip,
} from 'recharts';
import { lifecycle, compose, withState, branch, renderNothing, withProps } from 'recompose';
import { colsToRows } from '../webApp/src/common/utils';
import { graphql } from 'react-apollo';
import gql from 'graphql-tag';
import { connect } from 'react-redux';

const VizThumbNail = ({ isClient, coords }) => (
  <div style={{ height: '100%', width: '100%' }}>
    {isClient &&
      coords && (
        <ResponsiveContainer width="100%" height="100%">
          <ScatterChart margin={{ top: 0, right: 0, bottom: 0, left: 0 }}>
            <XAxis dataKey="dim1" type="number" hide />
            <YAxis dataKey="dim2" type="number" hide />
            <Scatter data={coords} fill="#8884d8" isAnimationActive={false} />
            <Tooltip cursor={{ strokeDasharray: '3 3' }} isAnimationActive={false} />
          </ScatterChart>
        </ResponsiveContainer>
      )}
  </div>
);

export const queryVizProbeKinds = gql`
  query VizProbeKinds($currentProjectRank: Int!) {
    whoami {
      id
      projectsByOwnerId(condition: { rank: $currentProjectRank }) {
        edges {
          node {
            id
            stepsByProjectId {
              edges {
                node {
                  id
                  exportsByStepId {
                    edges {
                      node {
                        id
                        kind
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
`;

const enhance = compose(
  withState('isClient', 'setIsClient', false),
  lifecycle({
    componentDidMount() {
      this.props.setIsClient(true);
    },
  }),
  connect(state => ({
    currentProjectRank: state.location.payload.projectRank,
  })),
  branch(props => props.currentProjectRank == null, renderNothing),
  graphql(queryVizProbeKinds),
  branch(({ data }) => data.loading, renderNothing),
  branch(({ data }) => data.whoami == null, renderNothing),
  branch(({ data }) => data.whoami.projectsByOwnerId.edges[0] == null, renderNothing),
  withProps(({ data }) => ({
    availableCoordses: data.whoami.projectsByOwnerId.edges[0].node.stepsByProjectId.edges
      .map(x => x.node.exportsByStepId.edges.map(y => ({ ...y.node })))
      .reduce((x, y) => x.concat(y || []), [])
      .filter(x => x.kind === 'sampleCoords'),
  })),
  branch(({ availableCoordses }) => availableCoordses.length === 0, renderNothing),
  withProps(({ availableCoordses }) => ({
    coordsId: availableCoordses[0].id,
  })),
  graphql(
    gql`
      query VizGetData($coordsId: UUID!) {
        exportById(id: $coordsId) {
          meta
          data
        }
      }
    `,
  ),
  branch(({ data }) => data.loading || data.exportById.data == null, renderNothing),
  withProps(({ data }) => ({
    coordsData: data.exportById.data,
  })),
  withProps(({ coordsData }) => ({
    coords: colsToRows({ dim1: coordsData[0].data, dim2: coordsData[1].data }),
  })),
);

export default enhance(VizThumbNail);
