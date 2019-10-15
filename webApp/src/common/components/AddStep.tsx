import { Avatar, List, ListItem, ListItemText, ListSubheader } from '@material-ui/core';
import gql from 'graphql-tag';
import React from 'react';
import { graphql } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, defaultProps, renderNothing, withProps } from 'recompose';

import actionCreators from '../redux/actionCreators';
import { guardEmpty } from '../utils';

// const Rating = ({ score, maxScore = 5 }) => {
//   const roundedScore = Math.round(score * 2) / 2;
//   const numFullStars = Math.floor(roundedScore);
//   const halfStar = roundedScore - Math.floor(roundedScore) > Number.EPSILON;
//   const numEmptyStars = Math.floor(maxScore - numFullStars - (halfStar ? 1 : 0));
//   return (
//     <div
//       style={{
//         flexGrow: 0,
//         flexShrink: 0,
//         flexBasis: 24 * maxScore,
//       }}
//     >
//       {[...Array(numFullStars)].map((e, i) => (
//         // eslint-disable-next-line react/no-array-index-key
//         <Star key={i} />
//       ))}
//       {halfStar && <StarHalf />}
//       {[...Array(numEmptyStars)].map((e, i) => (
//         // eslint-disable-next-line react/no-array-index-key
//         <StarBorder key={i} />
//       ))}
//     </div>
//   );
// };

const AddStep: React.FunctionComponent<any> = ({
  projectId,
  atStepRank,
  /**/
  rfrProjectStep,
  onClose,
  filteredGboxes,
  addStep,
}) => (
  <>
    <List subheader={<ListSubheader disableSticky>Individual steps</ListSubheader>}>
      {filteredGboxes.map((gbox) => (
        <ListItem
          button
          key={gbox.id}
          onClick={() => {
            addStep({
              gboxId: gbox.id,
              projectId,
              atStepRank,
              jumpToNewStep: true,
            });
            onClose();
          }}
        >
          <Avatar>{gbox.meta.title.slice(0, 2)}</Avatar>
          <ListItemText primary={gbox.meta.title} secondary={gbox.meta.subtitle} />
        </ListItem>
      ))}
    </List>
  </>
);

const enhance: any = compose(
  defaultProps({
    onClose: () => null,
  }),
  connect(
    null,
    {
      addStep: actionCreators.addStep,
      rfrProjectStep: actionCreators.rfrProjectStep,
    },
  ),
  graphql(
    gql`
      query GetGboxes {
        allGboxes {
          nodes {
            id
            meta
          }
        }
      }
    `,
    { options: { fetchPolicy: 'network-only' } },
  ),
  branch(({ data }) => data.loading, renderNothing),
  /* TODO: this is a hack, remove when upstream fixes the bug */
  guardEmpty('allGboxes'),
  withProps(({ data }) => ({
    gboxes: data.allGboxes.nodes,
  })),
  withProps(({ gboxes, filterText }) => ({
    filteredGboxes: gboxes.filter(
      (x) =>
        (x.meta.subtitle && x.meta.subtitle.toLowerCase().search(filterText.toLowerCase()) !== -1) ||
        (x.meta.title && x.meta.title.toLowerCase().search(filterText.toLowerCase()) !== -1),
    ),
  })),
);

export default enhance(AddStep);
