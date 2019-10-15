import gql from 'graphql-tag';
import { RFR_PROJECT_STEP, UPDATE_STEP_JUST_CHANGED_STATUS } from '../common/constants';
import { delay } from '../common/utils';

const PROBE_INTERVAL = 1000;

export default async (store, apolloClient) => {
  let prevStepProbes = null;

  // noinspection InfiniteLoopJS
  while (true) {
    const state = store.getState();

    if (state.location.type !== RFR_PROJECT_STEP) {
      await delay(PROBE_INTERVAL);
      continue;
    }
    // Check the UpdateAt for all steps in the current project.
    // If the server has a newer UpdateAt, update the data

    const currentProjectId = state.app.currentProjectId;
    const currentStepId = state.app.currentStepId;

    if (currentProjectId == null || currentStepId == null) {
      await delay(PROBE_INTERVAL);
      continue;
    }

    const resPollingWhetherNecessary = await apolloClient.query({
      query: gql`
        query PollingWhetherNecessary($currentProjectId: UUID!) {
          projectById(id: $currentProjectId) {
            id
            stepsByProjectId {
              nodes {
                id
                status
              }
            }
          }
        }
      `,
      variables: { currentProjectId },
    });

    const stepStatusesProject = resPollingWhetherNecessary.data.projectById;

    if (stepStatusesProject == null) {
      await delay(PROBE_INTERVAL);
      continue;
    }

    const stepStatuses = stepStatusesProject.stepsByProjectId.nodes.map((x) => x.status);

    if (
      stepStatuses.indexOf('INITIATED') === -1 &&
      stepStatuses.indexOf('INTERCEPTION_REQUESTED') === -1 &&
      stepStatuses.indexOf('RUNNING') === -1
    ) {
      await delay(PROBE_INTERVAL);
      continue;
    }

    const resPollingTick = await apolloClient.query({
      query: gql`
        query PollingGetUpdateAt($currentProjectId: UUID!) {
          projectById(id: $currentProjectId) {
            id
            stepsByProjectId {
              nodes {
                id
                updatedAt
              }
            }
          }
        }
      `,
      variables: { currentProjectId },
      fetchPolicy: 'network-only',
    });

    const stepProbes = resPollingTick.data.projectById.stepsByProjectId.nodes;

    if (prevStepProbes == null) {
      prevStepProbes = stepProbes;
      await delay(PROBE_INTERVAL);
      continue;
    }

    await Promise.all(
      stepProbes.map(async ({ id, status, updatedAt }) => {
        const prevStepProbe = prevStepProbes.find((x) => x.id === id);
        const prevUpdatedAt = prevStepProbe && prevStepProbe.updatedAt;

        const res1 = await apolloClient.query({
          query: gql`
            query PollingTick_getInfo($id: UUID!) {
              stepById(id: $id) {
                id
                status
              }
            }
          `,
          variables: {
            id,
          },
        });

        const prevStatus = res1.data.stepById.status;
        const stepId = res1.data.stepById.id;

        if (prevUpdatedAt === updatedAt) {
          return;
        }

        const res = await apolloClient.query({
          query: gql`
            query PollingUpdateStep($id: UUID!) {
              stepById(id: $id) {
                id
                status
                exportsByStepId {
                  nodes {
                    id
                    kind
                    meta
                    name
                  }
                }
                results
                errors
                state
              }
            }
          `,
          variables: { id },
          fetchPolicy: 'network-only',
        });

        const newStatus = res.data.stepById.status;

        if (currentStepId !== stepId && prevStatus !== 'DONE' && newStatus === 'DONE') {
          store.dispatch({
            type: UPDATE_STEP_JUST_CHANGED_STATUS,
            payload: { $set: { open: true, stepId, prevStatus, newStatus } },
          });
        }
      }),
    );

    prevStepProbes = stepProbes;
    await delay(PROBE_INTERVAL);
  }
};
