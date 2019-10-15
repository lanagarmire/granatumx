import {
  Avatar,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  List,
  ListItem,
  ListItemText,
  ListSubheader,
  Typography as T,
  withWidth,
} from '@material-ui/core';
import gql from 'graphql-tag';
import _ from 'lodash';
import React from 'react';
import { graphql } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, defaultProps, renderNothing, withProps } from 'recompose';

import actionCreators from '../redux/actionCreators';
import { withStateUpdater } from '../utils';
import Pipeline from './Pipeline';

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

const AddRecipe: React.FunctionComponent<any> = ({
  projectId,
  atStepRank,
  /**/
  numSteps,
  replaceOrAppendDialog,
  updateReplaceOrAppendDialog,
  removeAllSteps,
  addRecipeSteps,
  filteredRecipes,
  onClose,
  width,
}) => (
  <div>
    <List subheader={<ListSubheader disableSticky>Recipes</ListSubheader>}>
      {filteredRecipes.map((r) => (
        <ListItem
          key={r.id}
          button
          onClick={() => {
            if (numSteps > 0) {
              updateReplaceOrAppendDialog({ $set: { open: true, recipeId: r.id } });
            } else {
              addRecipeSteps({
                recipeId: r.id,
                projectId,
                atStepRank,
              });
              onClose();
            }
          }}
        >
          {(({ title, subtitle, steps }) => (
            <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'start' }}>
              <div style={{ display: 'flex', flexDirection: 'row', alignItems: 'center' }}>
                <Avatar>{title.slice(0, 2)}</Avatar>
                <ListItemText primary={title} secondary={subtitle} />
              </div>
              {width !== 'xs' && width !== 'sm' && <Pipeline steps={steps} />}
            </div>
          ))({
            title: r.meta.title,
            subtitle: r.meta.subtitle,
            steps: _(r.recipeGboxesByRecipeId.nodes)
              .orderBy('rank')
              .map((x) => ({ id: x.id, title: x.gboxByGboxId.meta.title }))
              .value(),
          })}
        </ListItem>
      ))}
    </List>
    <Dialog open={replaceOrAppendDialog.open}>
      <DialogContent>
        <T>Do you want to replace the current steps with these steps, or append them after current step?</T>
      </DialogContent>
      <DialogActions>
        <Button
          onClick={() => {
            updateReplaceOrAppendDialog({ open: { $set: false } });
          }}
        >
          Cancel
        </Button>
        <Button
          onClick={() => {
            addRecipeSteps({
              recipeId: replaceOrAppendDialog.recipeId,
              projectId,
              atStepRank,
            });
            updateReplaceOrAppendDialog({ open: { $set: false } });
            onClose();
          }}
        >
          Append
        </Button>
        <Button
          color="primary"
          onClick={() => {
            removeAllSteps({
              projectId,
              cb: () => {
                addRecipeSteps({
                  recipeId: replaceOrAppendDialog.recipeId,
                  projectId,
                  atStepRank,
                });
              },
            });
            updateReplaceOrAppendDialog({ open: { $set: false } });
            onClose();
          }}
        >
          Replace
        </Button>
      </DialogActions>
    </Dialog>
  </div>
);

const enhance: any = compose(
  defaultProps({
    onClose: () => null,
  }),
  connect(
    null,
    {
      addRecipeSteps: actionCreators.addRecipeSteps,
      removeAllSteps: actionCreators.removeAllSteps,
    },
  ),
  graphql(
    gql`
      query GetGboxes($projectId: UUID!) {
        projectById(id: $projectId) {
          id
          stepsByProjectId {
            totalCount
          }
        }
        allRecipes {
          nodes {
            id
            meta
            recipeGboxesByRecipeId {
              nodes {
                id
                rank
                gboxId
                gboxByGboxId {
                  id
                  meta
                }
              }
            }
          }
        }
      }
    `,
    { options: { fetchPolicy: 'network-only' } },
  ),
  branch(({ data }) => data.loading, renderNothing),
  withProps(({ data }) => ({
    recipes: _(data.allRecipes.nodes)
      .orderBy(['meta.featured', 'meta.title'])
      .value(),
    numSteps: data.projectById.stepsByProjectId.totalCount,
  })),
  withProps(({ recipes, filterText }) => ({
    filteredRecipes: recipes.filter(
      (x) =>
        x.meta.subtitle.toLowerCase().search(filterText.toLowerCase()) !== -1 ||
        x.meta.title.toLowerCase().search(filterText.toLowerCase()) !== -1 ||
        x.recipeGboxesByRecipeId.nodes.some(
          (y) => y.gboxByGboxId.meta.title.toLowerCase().search(filterText.toLowerCase()) !== -1,
        ),
    ),
  })),
  withStateUpdater('replaceOrAppendDialog', { open: false, recipeId: undefined }),
  withWidth(),
);

export default enhance(AddRecipe);
