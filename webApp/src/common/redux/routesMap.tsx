import gql from 'graphql-tag';
import update from 'immutability-helper';
import _ from 'lodash';
import { NOT_FOUND, redirect } from 'redux-first-router';
import { queryBackend, trace } from '../utils';

import {
  ADD_RECIPE_STEPS,
  ADD_STEP,
  CLEAR_ERRORS_FOR_STEP,
  CREATE_PROJECT,
  DELETE_PROJECT,
  DYNAMIC_DISPATCH_ACTION,
  GLOBAL_ERROR,
  INTERCEPT_STEP,
  MODIFY_PROJECT,
  NEXT_STEP,
  PREV_STEP,
  QUERY_BACKEND,
  REMOVE_ALL_STEPS,
  REMOVE_STEP,
  REORDER_STEP,
  RESET_STEP,
  RESET_STEP_RECURSIVELY,
  RFR_CREATE_PROJECT,
  RFR_HOME,
  RFR_JUMP_TO_STEP_BY_STEP_RANK,
  RFR_LOGOUT,
  RFR_NO_STEPS,
  RFR_NO_STEPS_THUNK_END,
  RFR_NOT_FOUND,
  RFR_PROJECT_REPORT,
  RFR_PROJECT_ROOT,
  RFR_PROJECT_STEP,
  RFR_PROJECT_STEP_THUNK_END,
  SAVE_CURRENT_PIPELINE_AS_RECIPE,
  SUBMIT_STEP,
  UPDATE_CURRENT_STEP,
  UPDATE_WELCOME_DIALOG,
  UPLOAD_FILE,
} from '../constants';

const queryFindLastStep = gql`
  query FindLastStep($projectRank: Int) {
    whoami {
      id
      projectsByOwnerId(condition: { rank: $projectRank }, orderBy: RANK_DESC, first: 1) {
        nodes {
          id
          rank
          stepsByProjectId(orderBy: RANK_ASC) {
            nodes {
              id
              rank
              status
            }
          }
        }
      }
    }
  }
`;

const download = (filename, text) => {
  if (typeof document === 'object') {
    const element = document.createElement('a');
    element.setAttribute('href', `data:text/plain;charset=utf-8,${encodeURIComponent(text)}`);
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
  } else if (__DEV__) {
    throw new Error('Trying to download on the *server-side*. Rejected.');
  }
};

export default {
  // TODO: configure redirection
  [RFR_HOME]: {
    path: '/',
    thunk: async (dispatch, getState, { extra: { apolloClient } }) => {
      try {
        const res = await apolloClient.query({
          query: gql`
            query RFR_HOME {
              whoami {
                id
                state
              }
            }
          `,
        });

        const currentStepId = res.data.whoami.state && res.data.whoami.state.currentStepId;

        if (currentStepId != null) {
          dispatch(redirect({ type: RFR_PROJECT_STEP, payload: { stepId: currentStepId } }));
          return;
        }

        const resGetFirstProject = await apolloClient.query({
          query: gql`
            query RFR_HOME_GetFirstProject {
              whoami {
                id
                projectsByOwnerId(orderBy: NAME_ASC, first: 1) {
                  nodes {
                    id
                  }
                }
              }
            }
          `,
        });

        const firstProject = resGetFirstProject.data.whoami.projectsByOwnerId.nodes[0];

        if (firstProject != null) {
          dispatch(redirect({ type: RFR_PROJECT_ROOT, payload: { projectId: firstProject.id } }));
          return;
        }

        // If you are here, you don't even have a project
        dispatch(
          redirect({
            type: CREATE_PROJECT,
            payload: { name: 'New sandbox project', description: 'Add description here.', welcome: true },
          }),
        );
      } catch (e) {
        if (__DEV__) {
          console.error(e);
        }
        dispatch(redirect({ type: RFR_NOT_FOUND }));
      }
    },
  },

  [RFR_CREATE_PROJECT]: {
    path: '/create-project',
  },

  [RFR_NOT_FOUND]: {
    thunk: async (dispatch, getState, { action: { payload }, extra: { apolloClient } }) => {
      dispatch({ type: UPDATE_CURRENT_STEP, payload: { stepId: null } });
      dispatch(redirect({ type: NOT_FOUND }));
    },
  },

  [CREATE_PROJECT]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { name, description, welcome },
        },
        extra: { apolloClient },
      },
    ) => {
      const res = await apolloClient.mutate({
        mutation: gql`
          mutation CREATE_PROJECT($name: String, $description: String) {
            createProject(input: { name: $name, description: $description }) {
              createProjectReturning {
                id
              }
              query {
                whoami {
                  id
                  projectsByOwnerId {
                    nodes {
                      id
                      rank
                      name
                      description
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { name, description },
      });

      const newProjectId = res.data.createProject.createProjectReturning.id;

      dispatch(redirect({ type: RFR_PROJECT_ROOT, payload: { projectId: newProjectId, welcome } }));
    },
  },

  [RFR_JUMP_TO_STEP_BY_STEP_RANK]: {
    path: '/project/:projectId/step-rank/:stepRank',
    toPath: (value, key) => (key === 'stepRank' ? (+value + 1).toString() : value),
    fromPath: (pathSegment, key) => (key === 'stepRank' ? +pathSegment - 1 : pathSegment),
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { projectId, stepRank },
        },
        extra: { apolloClient },
      },
    ) => {
      try {
        const res = await apolloClient.query({
          query: gql`
            query RFR_JUMP_TO_STEP_BY_STEP_RANK($projectId: UUID!, $stepRank: Int!) {
              projectById(id: $projectId) {
                id
                stepsByProjectId(condition: { rank: $stepRank }) {
                  nodes {
                    id
                  }
                }
              }
            }
          `,
          variables: { projectId, stepRank },
          // fetchPolicy: 'network-only',
        });

        const project = res.data.projectById;

        if (project == null) {
          throw new Error('Project not found.');
        }

        const step = project.stepsByProjectId.nodes[0];

        if (step == null) {
          throw new Error('Step not found');
        }

        dispatch({ type: RFR_PROJECT_STEP, payload: { stepId: step.id } });
      } catch (e) {
        if (__DEV__) {
          console.error(e);
        }
        dispatch(redirect({ type: RFR_NOT_FOUND }));
      }
    },
  },

  [RFR_PROJECT_ROOT]: {
    // path: '/project/:projectId',
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { projectId, welcome },
        },
        extra: { apolloClient },
      },
    ) => {
      try {
        const res = await apolloClient.query({
          query: gql`
            query RFR_PROJECT_ROOT($projectId: UUID!) {
              projectById(id: $projectId) {
                id
                state
              }
            }
          `,
          variables: { projectId },
        });

        const currentStepId = res.data.projectById.state && res.data.projectById.state.currentStepId;

        if (currentStepId != null) {
          dispatch(redirect({ type: RFR_PROJECT_STEP, payload: { stepId: currentStepId } }));
          return;
        }

        const resGetSteps = await apolloClient.query({
          query: gql`
            query RFR_PROJECT_ROOT_GetSteps($projectId: UUID!) {
              projectById(id: $projectId) {
                id
                stepsByProjectId {
                  totalCount
                }
              }
            }
          `,
          variables: { projectId },
        });

        const numSteps = resGetSteps.data.projectById.stepsByProjectId.totalCount;

        if (welcome) {
          dispatch({ type: UPDATE_WELCOME_DIALOG, payload: { open: { $set: true } } });
        }

        if (numSteps !== 0) {
          dispatch(redirect({ type: RFR_JUMP_TO_STEP_BY_STEP_RANK, payload: { projectId, stepRank: 0 } }));
          return;
        }

        dispatch(redirect({ type: RFR_NO_STEPS, payload: { projectId, stepRank: 0 } }));
      } catch (e) {
        if (__DEV__) {
          if (__DEV__) {
            console.error(e);
          }
        }
        dispatch(redirect({ type: RFR_NOT_FOUND }));
      }
    },
  },

  [RFR_PROJECT_REPORT]: {
    path: '/project/:projectId/report',
  },

  [RFR_NO_STEPS]: {
    path: '/project/:projectId/no-steps',
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { projectId },
        },
        extra: { apolloClient },
      },
    ) => {
      dispatch({ type: UPDATE_CURRENT_STEP, payload: { stepId: null } });

      const res = await apolloClient.query({
        query: gql`
          query RFR_NO_STEPS($projectId: UUID!) {
            projectById(id: $projectId) {
              id
              rank
              name
            }
          }
        `,
        variables: { projectId },
      });

      const project = res.data.projectById;

      if (project == null) {
        dispatch(redirect({ type: RFR_NOT_FOUND }));
        return;
      }

      dispatch({ type: RFR_NO_STEPS_THUNK_END, payload: { project } });
    },
  },

  [RFR_PROJECT_STEP]: {
    path: '/step/:stepId',
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      try {
        const res = await apolloClient.query({
          query: gql`
            query RFR_PROJECT_STEP($stepId: UUID!) {
              whoami {
                id
                profile {
                  displayName
                }
              }
              stepById(id: $stepId) {
                id
                rank
                projectByProjectId {
                  id
                  rank
                  name
                }
                gboxByGbox {
                  id
                  meta
                }
              }
            }
          `,
          variables: { stepId },
          // fetchPolicy: 'network-only',
        });

        const whoami = res.data.whoami;

        if (whoami == null) {
          throw new Error('You do not have a user id.');
        }

        const step = res.data.stepById;

        if (step == null) {
          throw new Error('Step not found.');
        }

        const project = step.projectByProjectId;

        if (project == null) {
          throw new Error('This should be impossible but: a step was found but there it has no owning project');
        }

        const gbox = step.gboxByGbox;

        if (gbox == null) {
          throw new Error('This should be impossible but: a step has not gbox');
        }

        // updating the current step *asynchronously*
        dispatch({ type: UPDATE_CURRENT_STEP, payload: { stepId } });

        dispatch({
          type: RFR_PROJECT_STEP_THUNK_END,
          payload: { whoami, project, step },
        });
      } catch (e) {
        if (__DEV__) {
          console.error(e);
        }
        dispatch(redirect({ type: RFR_NOT_FOUND }));
      }
    },
  },

  [PREV_STEP]: {
    thunk: async (dispatch, getState, { extra: { apolloClient } }) => {
      const currentStepId = getState().app.currentStepId;
      if (currentStepId == null) return;

      const res = await apolloClient.query({
        query: gql`
          query PREV_STEP($stepId: UUID!) {
            stepById(id: $stepId) {
              id
              rank
              projectByProjectId {
                id
                stepsByProjectId {
                  nodes {
                    id
                    rank
                  }
                }
              }
            }
          }
        `,
        variables: { stepId: currentStepId },
      });

      const currentStep = res.data.stepById;
      const currentProject = currentStep.projectByProjectId;
      const steps = currentProject.stepsByProjectId.nodes;

      const nextStep = steps.find((s) => s.rank === currentStep.rank - 1);
      if (nextStep != null) {
        dispatch({ type: RFR_PROJECT_STEP, payload: { stepId: nextStep.id } });
      }
    },
  },

  [NEXT_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { reportAtTheEnd = false },
        },
        extra: { apolloClient },
      },
    ) => {
      const currentStepId = getState().app.currentStepId;
      if (currentStepId == null) return;

      const res = await apolloClient.query({
        query: gql`
          query NEXT_STEP($stepId: UUID!) {
            stepById(id: $stepId) {
              id
              rank
              projectByProjectId {
                id
                stepsByProjectId {
                  nodes {
                    id
                    rank
                  }
                }
              }
            }
          }
        `,
        variables: { stepId: currentStepId },
      });

      const currentStep = res.data.stepById;
      const currentProject = currentStep.projectByProjectId;
      const steps = currentProject.stepsByProjectId.nodes;

      const nextStep = steps.find((s) => s.rank === currentStep.rank + 1);
      console.log('reportAtTheEnd =', reportAtTheEnd);
      if (reportAtTheEnd && nextStep == null) {
        dispatch({ type: RFR_PROJECT_REPORT, payload: { projectId: currentProject.id } });
      } else {
        dispatch({ type: RFR_PROJECT_STEP, payload: { stepId: nextStep.id } });
      }
    },
  },

  [RFR_LOGOUT]: {
    path: '/logout',
  },

  [UPDATE_CURRENT_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      // TODO: upate currentStepId for the current *Project* too
      const resOldStates = await apolloClient.query({
        query: gql`
          query UPDATE_CURRENT_STEP_OldStates {
            whoami {
              id
              state
            }
          }
        `,
      });

      const oldWhoamiState = resOldStates.data.whoami.state;
      const newWhoamiState = update(oldWhoamiState || {}, { currentStepId: { $set: stepId } });

      await apolloClient.mutate({
        mutation: gql`
          mutation UPDATE_CURRENT_STEP_UpdateState($newWhoamiState: JSON!) {
            updateWhoamiState(input: { newState: $newWhoamiState }) {
              query {
                whoami {
                  id
                  state
                }
              }
            }
          }
        `,
        variables: { newWhoamiState },
      });
    },
  },

  [MODIFY_PROJECT]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { id, name, description },
        },
        extra: { apolloClient },
      },
    ) => {
      await apolloClient.mutate({
        mutation: gql`
          mutation MODIFY_PROJECT($id: UUID!, $name: String, $description: String) {
            modifyProject(id: $id, name: $name, description: $description) {
              query {
                projectById(id: $id) {
                  id
                  name
                  description
                }
              }
            }
          }
        `,
        variables: { id, name, description },
      });
    },
  },

  [DELETE_PROJECT]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { projectId },
        },
        extra: { apolloClient },
      },
    ) => {
      if (!confirm('Are you sure?\n\nDeleting this project will REMOVE all data in it.')) {
        return;
      }

      await apolloClient.mutate({
        mutation: gql`
          mutation DELETE_PROJECT($projectId: UUID!) {
            deleteProject(input: { projectId: $projectId }) {
              query {
                whoami {
                  id
                  projectsByOwnerId {
                    nodes {
                      id
                      rank
                      name
                      description
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { projectId },
      });

      if (typeof window === 'object') {
        window.location.href = '/';
      }
    },
  },

  [UPLOAD_FILE]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { file, onProgress, onLoad, onError, onResponse },
        },
        extra: { fetch },
      },
    ) => {
      const fd = new FormData();
      fd.append('file', file);

      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/file-upload', true);
      xhr.upload.onprogress = onProgress;
      xhr.upload.onload = onLoad;
      xhr.onload = () => {
        onResponse(xhr.response);
      };
      xhr.upload.onerror = onError;
      xhr.send(fd);
    },
  },

  [QUERY_BACKEND]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { cb, ...payload },
        },
        extra: { fetch },
      },
    ) => {
      const results = await queryBackend(fetch, payload);
      if (typeof cb === 'function') {
        cb(results);
      }
    },
  },

  [RESET_STEP_RECURSIVELY]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      try {
        await apolloClient.mutate({
          mutation: gql`
            mutation RESET_STEP_RECURSIVELY($stepId: UUID!) {
              resetStepRecursively(input: { stepId: $stepId }) {
                query {
                  stepById(id: $stepId) {
                    id
                    projectByProjectId {
                      id
                      stepsByProjectId {
                        nodes {
                          id
                          status
                          args
                          importsByStepId {
                            nodes {
                              id
                              stepId
                              exportId
                              injectInto
                            }
                          }
                          exportsByStepId {
                            nodes {
                              id
                              stepId
                              extractFrom
                              kind
                              meta
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          `,
          variables: { stepId },
        });
      } catch (e) {
        if (__DEV__) {
          console.error(e);
        }
        dispatch(redirect({ type: RFR_NOT_FOUND }));
      }
    },
  },

  [RESET_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      try {
        if (!confirm('Are you sure?\n\nResetting this step will remove all the data generated by it.')) {
          return;
        }

        await apolloClient.mutate({
          mutation: gql`
            mutation RESET_STEP($stepId: UUID!) {
              resetStep(input: { stepId: $stepId }) {
                query {
                  stepById(id: $stepId) {
                    id
                    status
                    args
                    importsByStepId {
                      nodes {
                        id
                        stepId
                        exportId
                        injectInto
                      }
                    }
                    exportsByStepId {
                      nodes {
                        id
                        stepId
                        extractFrom
                        kind
                        meta
                      }
                    }
                  }
                }
              }
            }
          `,
          variables: { stepId },
        });
      } catch (e) {
        dispatch({
          type: GLOBAL_ERROR,
          payload: {
            title: 'Cannot reset',
            content: e.message,
            actions: [
              {
                label: 'Force reset',
                description: 'Reset all the steps that depend on this step, and then reset this step',
                action: { type: RESET_STEP_RECURSIVELY, payload: { stepId } },
              },
            ],
          },
        });
        throw e;
      }
    },
  },

  [DYNAMIC_DISPATCH_ACTION]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { action },
        },
      },
    ) => {
      dispatch(action);
    },
  },

  // When user clicks the ">" button of a step
  [SUBMIT_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { id, gbox, args, imports, uploadedFiles, exports, results },
        },
        extra: { apolloClient, fetch },
      },
    ) => {
      // const { location: { payload: { stepRank, projectRank } } } = getState();
      // save the step before doing anything else
      // TODO: change all clientMutationId into the things expected to be changed

      await queryBackend(fetch, {
        audience: '__granatum',
        endpoint: 'saveStep',
        stepId: id,
        args,
        imports,
        uploadedFiles,
        exports,
        results,
      });

      await apolloClient.query({
        query: gql`
          query SUBMIT_STEP_save($id: UUID!) {
            stepById(id: $id) {
              id
              status
              args
              results
              errors
              state
              importsByStepId {
                nodes {
                  id
                  stepId
                  exportId
                  injectInto
                }
              }
              exportsByStepId {
                nodes {
                  id
                  stepId
                  extractFrom
                  kind
                  meta
                }
              }
            }
          }
        `,
        variables: { id },
        fetchPolicy: 'network-only',
      });

      // If gbox is null it means it's a client-only module (like file upload)
      if (gbox == null) {
        console.log('This is a client-only module. No back-end is invoked.');

        await apolloClient.mutate({
          mutation: gql`
            mutation SUBMIT_STEP_markDone($input: MarkDoneStepInput!, $id: UUID!) {
              markDoneStep(input: $input) {
                query {
                  stepById(id: $id) {
                    id
                    status
                    exportsByStepId {
                      nodes {
                        id
                      }
                    }
                  }
                }
              }
            }
          `,
          variables: { input: { id }, id },
        });

        return;
      }

      // otherwise, change the status to "initiated" and wait for the task scheduler to pick it up
      await apolloClient.mutate({
        mutation: gql`
          mutation SUBMIT_STEP_initiate($input: InitiateStepInput!, $id: UUID!) {
            initiateStep(input: $input) {
              query {
                stepById(id: $id) {
                  id
                  status # expected to be "initiated" and shortly "running"
                  # Note that exports are updated immediately without waiting for the step to finish
                  # its back-end execution. This is so that users can start working on the subsequent steps
                  # using the un-populated exports as dependencies immediately
                  exportsByStepId {
                    nodes {
                      id
                      kind
                      meta
                      extractFrom
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { input: { id }, id },
      });

      // TODO: Here, we should start a polling loop getting the newest status of this step, a spinner
      // TODO:     should be shown as well.
    },
  },

  [ADD_RECIPE_STEPS]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { recipeId, projectId, atStepRank },
        },
        extra: { apolloClient, fetch },
      },
    ) => {
      const resGetStepsFromRecipe = await apolloClient.query({
        query: gql`
          query ADD_RECIPE_STEPS_GetStepsFromRecipe($recipeId: UUID!) {
            recipeById(id: $recipeId) {
              id
              recipeGboxesByRecipeId {
                nodes {
                  id
                  rank
                  gboxId
                  initialState
                }
              }
            }
          }
        `,
        variables: { recipeId },
      });

      const gboxes = _(resGetStepsFromRecipe.data.recipeById.recipeGboxesByRecipeId.nodes)
        .orderBy('rank')
        .map((x) => ({
          gboxId: x.gboxId,
          state: x.initialState,
        }))
        .value();

      await apolloClient.mutate({
        mutation: gql`
          mutation ADD_RECIPE_STEPS($gboxes: [AddStepsGboxInput]!, $projectId: UUID!, $atStepRank: Int!) {
            addSteps(input: { atStepRank: $atStepRank, gboxes: $gboxes, projectId: $projectId }) {
              query {
                projectById(id: $projectId) {
                  id
                  stepsByProjectId {
                    nodes {
                      id
                      rank
                      status
                      gboxByGbox {
                        id
                        meta
                      }
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { gboxes, projectId, atStepRank },
      });
    },
  },

  [ADD_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { gboxId, projectId, atStepRank, jumpToNewStep },
        },
        extra: { apolloClient, fetch },
      },
    ) => {
      const res = await apolloClient.mutate({
        mutation: gql`
          mutation ADD_STEP($gboxId: String!, $projectId: UUID!, $atStepRank: Int!) {
            addStep(input: { atStepRank: $atStepRank, gboxId: $gboxId, projectId: $projectId }) {
              addStepReturning {
                newStepId
              }
              query {
                projectById(id: $projectId) {
                  id
                  stepsByProjectId {
                    nodes {
                      id
                      rank
                      status
                      gboxByGbox {
                        id
                        meta
                      }
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { gboxId, projectId, atStepRank },
      });

      const newStepId = res.data.addStep.addStepReturning.newStepId;

      if (jumpToNewStep) {
        dispatch({ type: RFR_PROJECT_STEP, payload: { stepId: newStepId } });
      }
    },
  },

  [REORDER_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId, toRank, cb },
        },
        extra: { apolloClient, fetch },
      },
    ) => {
      const resGetCurrent = await apolloClient.query({
        query: gql`
          query REORDER_STEP_GetCurrent($stepId: UUID!) {
            stepById(id: $stepId) {
              id
              projectByProjectId {
                id
                stepsByProjectId {
                  nodes {
                    id
                    rank
                    status
                    gboxByGbox {
                      id
                      meta
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { stepId },
      });

      const affectedProjectId = resGetCurrent.data.stepById.projectByProjectId.id;

      try {
        await apolloClient.mutate({
          mutation: gql`
            mutation REORDER_STEP($stepId: UUID!, $toRank: Int!, $affectedProjectId: UUID!) {
              reorderStep(input: { id: $stepId, toRank: $toRank }) {
                query {
                  projectById(id: $affectedProjectId) {
                    id
                    stepsByProjectId {
                      nodes {
                        id
                        rank
                        status
                        gboxByGbox {
                          id
                          meta
                        }
                      }
                    }
                  }
                }
              }
            }
          `,
          // optimisticResponse: {
          //   __typename: 'Mutation',
          //   removeStep: {
          //     __typename: 'RemoveStepPayload',
          //     query: {
          //       __typename: 'Query',
          //       stepsByProjectId: {
          //         ...update(resGetCurrent.data.stepById.stepsByProjectId, {
          //           nodes: (nodes) => {
          //             const rankToRemove = nodes.find((x) => x.id === stepId).rank;
          //             return nodes
          //               .filter((x) => x.id !== stepId)
          //               .map((x) => (x.rank > rankToRemove ? update(x, { rank: (x) => x - 1 }) : x));
          //           },
          //         }),
          //       },
          //     },
          //   },
          // },
          variables: { stepId, toRank, affectedProjectId },
        });
      } catch (e) {
        dispatch({
          type: GLOBAL_ERROR,
          payload: {
            title: 'Cannot reorder',
            content: e.message,
          },
        });
        throw e;
      }

      if (typeof cb === 'function') {
        cb();
      }
    },
  },

  [REMOVE_ALL_STEPS]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { projectId, cb },
        },
        extra: { apolloClient },
      },
    ) => {
      // TODO: if any step is running, throw an error
      await apolloClient.mutate({
        mutation: gql`
          mutation REMOVE_ALL_STEPS($projectId: UUID!) {
            clearStepsInProject(input: { projectId: $projectId }) {
              query {
                projectById(id: $projectId) {
                  id
                  stepsByProjectId {
                    nodes {
                      id
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { projectId },
      });

      if (typeof cb === 'function') {
        cb();
      }
    },
  },

  [REMOVE_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId, doRedirect },
        },
        extra: { apolloClient, fetch },
      },
    ) => {
      // const checkRevDepsRes = await apolloClient.query({
      //   query: gql`
      //     query checkRevDeps($stepId: UUID!) {
      //       stepById {
      //         id
      //         status
      //       }
      //       getStepRevDeps(id: $stepId) {
      //         nodes {
      //           id
      //           rank
      //           gboxByGbox {
      //             id
      //             meta
      //           }
      //         }
      //       }
      //     }
      //   `,
      //   variables: { stepId },
      //   fetchPolicy: 'network-only',
      // });
      //
      // const revDeps = checkRevDepsRes.data.getStepRevDeps.nodes.map((x) => ({
      //   rank: x.rank,
      //   title: x.gboxByGbox.meta.title,
      // }));
      //
      // if (revDeps.length > 0) {
      //   dispatch({
      //     type: GLOBAL_ERROR,
      //     payload: {
      //       title: 'Cannot remove step',
      //       content: [
      //         'The step you are trying to remove is depended upon by the following step(s).',
      //         'You have to remove them first:\n',
      //         '\n',
      //         ...revDeps.map((x) => `  * Step ${x.rank + 1}: **${x.title}**\n`),
      //       ].join(''),
      //     },
      //   });
      //   return;
      // }

      // This is for optimistic response and updating the parent project
      const resGetCurrent = await apolloClient.query({
        query: gql`
          query REMOVE_STEP_GetCurrent($stepId: UUID!) {
            stepById(id: $stepId) {
              id
              projectByProjectId {
                id
                stepsByProjectId {
                  nodes {
                    id
                    rank
                    status
                    gboxByGbox {
                      id
                      meta
                    }
                  }
                }
              }
            }
          }
        `,
        variables: { stepId },
      });

      const affectedProjectId = resGetCurrent.data.stepById.projectByProjectId.id;
      const currentStatus = resGetCurrent.data.stepById.projectByProjectId.stepsByProjectId.nodes.find(
        (x) => x.id === stepId,
      ).status;

      if (currentStatus !== 'IDLE') {
        if (!confirm('Are you sure?\n\nRemoving this step will remove all the data generated by it.')) {
          return;
        }
      }

      try {
        await apolloClient.mutate({
          mutation: gql`
            mutation REMOVE_STEP($stepId: UUID!, $affectedProjectId: UUID!) {
              removeStep(input: { stepId: $stepId }) {
                query {
                  projectById(id: $affectedProjectId) {
                    id
                    stepsByProjectId {
                      nodes {
                        id
                        rank
                        status
                        gboxByGbox {
                          id
                          meta
                        }
                      }
                    }
                  }
                }
              }
            }
          `,
          optimisticResponse: {
            __typename: 'Mutation',
            removeStep: {
              __typename: 'RemoveStepPayload',
              query: {
                __typename: 'Query',
                projectById: {
                  ...update(resGetCurrent.data.stepById.projectByProjectId, {
                    stepsByProjectId: {
                      nodes: (nodes) => {
                        const rankToRemove = nodes.find((x) => x.id === stepId).rank;
                        return nodes
                          .filter((x) => x.id !== stepId)
                          .map((x) => (x.rank > rankToRemove ? update(x, { rank: (y) => y - 1 }) : x));
                      },
                    },
                  }),
                },
              },
            },
          },
          variables: { stepId, affectedProjectId },
        });
      } catch (e) {
        dispatch({
          type: GLOBAL_ERROR,
          payload: {
            title: 'Cannot remove this step',
            content: e.message,
          },
        });
        throw e;
      }

      if (doRedirect) {
        if (doRedirect.currentStepId == null || stepId === doRedirect.currentStepId) {
          dispatch({ type: RFR_NO_STEPS, payload: { projectId: doRedirect.currentProjectId } });
        }
      }
    },
  },

  [CLEAR_ERRORS_FOR_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      await apolloClient.mutate({
        mutation: gql`
          mutation CLEAR_ERRORS_FOR_STEP($input: ClearErrorsForStepInput!, $stepId: UUID!) {
            clearErrorsForStep(input: $input) {
              query {
                stepById(id: $stepId) {
                  id
                  errors
                }
              }
            }
          }
        `,
        variables: { input: { stepId }, stepId },
      });
    },
  },

  [SAVE_CURRENT_PIPELINE_AS_RECIPE]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { meta },
        },
        extra: { apolloClient },
      },
    ) => {
      const { componentStates, currentProjectId } = getState().app;

      const res = await apolloClient.query({
        query: gql`
          query SAVE_CURRENT_PIPELINE_AS_RECIPE_GetCurrentPipeline($currentProjectId: UUID!) {
            projectById(id: $currentProjectId) {
              id
              stepsByProjectId {
                nodes {
                  id
                  rank
                  state
                  gbox
                }
              }
            }
          }
        `,
        variables: { currentProjectId },
      });

      const gboxes = res.data.projectById.stepsByProjectId.nodes.map((s) => ({
        gboxId: s.gbox,
        rank: s.rank,
        // initialState: componentStates[s.id] != null ? componentStates[s.id] : s.state,
        initialState: null,
      }));

      await apolloClient.mutate({
        mutation: gql`
          mutation SAVE_CURRENT_PIPELINE_AS_RECIPE($meta: JSON!, $gboxes: [SaveCurrentPipelineAsRecipeGboxInput]!) {
            saveCurrentPipelineAsRecipe(input: { meta: $meta, gboxes: $gboxes }) {
              query {
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
            }
          }
        `,
        variables: { meta, gboxes },
      });
    },
  },

  [INTERCEPT_STEP]: {
    thunk: async (
      dispatch,
      getState,
      {
        action: {
          payload: { stepId },
        },
        extra: { apolloClient },
      },
    ) => {
      await apolloClient.mutate({
        mutation: gql`
          mutation INTERCEPT_STEP($stepId: UUID!) {
            interceptStep(input: { stepId: $stepId }) {
              query {
                stepById(id: $stepId) {
                  id
                  status
                }
              }
            }
          }
        `,
        variables: { stepId },
      });
    },
  },

  MUH_THUNK: {
    thunk: async (dispatch, getState, { action, extra }) => {
      await new Promise((resolve, reject) => {
        setTimeout(() => {
          resolve();
        }, 3000);
      });
      alert('Hello!');
    },
  },
};
