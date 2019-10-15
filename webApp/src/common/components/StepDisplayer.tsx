import { Typography as T } from '@material-ui/core';
import gql from 'graphql-tag';
import React from 'react';
import { graphql } from 'react-apollo';
import { connect } from 'react-redux';
import { branch, compose, lifecycle, renderNothing, withProps, withState } from 'recompose';
import actionCreators from '../redux/actionCreators';
import { guardEmpty } from '../utils';
import FrontEndFromSpecs from './FrontEndFromSpecs';
import InteractiveOutlierRemoval from './InteractiveOutlierRemoval';

const getCustomFrontend = (id) => {
  const customFrontends = {
    InteractiveOutlierRemoval,
  };

  const comp = customFrontends[id];

  if (comp == null) {
    return () => <T variant="caption">There is no custom frontend with id: {id}</T>;
  }

  return comp;
};

class ErrorBoundary extends React.Component<any, any> {
  constructor(props) {
    super(props);
    this.state = { hasError: false };
  }

  componentDidCatch(error, info) {
    // Display fallback UI
    this.setState({ hasError: true });
  }

  render() {
    if (this.state.hasError) {
      // You can render any custom fallback UI
      return <h1>Something went wrong. Please refresh your browser.</h1>;
    }
    return this.props.children;
  }
}

const StepDisplayer = ({ gFuncs, gboxId, meta, frontend }) => (
  <>
    {/* <JSONTree data={{ ...gFuncs, gboxId, meta, frontend }} /> */}
    {typeof frontend === 'object' ? (
      <ErrorBoundary>
        <FrontEndFromSpecs {...{ ...gFuncs, gboxId, meta, frontend }} />
      </ErrorBoundary>
    ) : typeof frontend === 'string' ? (
      React.createElement(getCustomFrontend(frontend), {
        ...gFuncs,
        gboxId,
        meta,
        frontend,
      })
    ) : (
      <T variant="caption">Unknown type of frontend: {typeof frontend}</T>
    )}
  </>
);

const enhance = compose<{ gFuncs; gboxId; meta; frontend }, any>(
  branch<{ currentStepId: string }>((props) => props.currentStepId == null, renderNothing),
  branch<{ currentProjectId: string }>((props) => props.currentProjectId == null, renderNothing),
  withState('isClient', 'setIsClient', false),
  lifecycle<{ setIsClient: (newIsClient: boolean) => null }, any, any>({
    componentDidMount() {
      this.props.setIsClient(true);
    },
  }),
  connect(
    null,
    {
      submitStep: actionCreators.submitStep,
      resetStep: actionCreators.resetStep,
      clearErrors: actionCreators.clearErrors,
      globalError: actionCreators.globalError,
      rfrProjectStep: actionCreators.rfrProjectStep,
      uploadFile: actionCreators.uploadFile,
      queryBackend: actionCreators.queryBackend,
      interceptStep: actionCreators.interceptStep,
      rfrProjectReport: actionCreators.rfrProjectReport,
      nextStep: actionCreators.nextStep,
    },
  ),
  graphql(
    gql`
      query ModuleChooser($currentProjectId: UUID!, $currentStepId: UUID!) {
        projectById(id: $currentProjectId) {
          id
          rank
          stepsByProjectId {
            nodes {
              id
              rank
              gboxByGbox {
                id
                meta
                frontend
              }
              exportsByStepId {
                nodes {
                  id
                  extractFrom
                  kind
                  meta
                }
              }
            }
          }
        }
        stepById(id: $currentStepId) {
          id
          rank
          status
          args
          results
          errors
          state
          gboxByGbox {
            id
            meta
            frontend
          }
          importsByStepId {
            nodes {
              id
              injectInto
              exportByExportId {
                id
                extractFrom
                stepByStepId {
                  id
                  rank
                  gboxByGbox {
                    id
                    meta
                  }
                }
              }
            }
          }
          uploadedFilesByStepId {
            nodes {
              id
              meta
              injectInto
            }
          }
          exportsByStepId {
            nodes {
              id
              extractFrom
              kind
              meta
            }
          }
        }
      }
    `,
  ),
  branch(({ data }) => data.loading, renderNothing),
  /* TODO: this is a hack, remove when upstream fixes the bug */
  guardEmpty('stepById'),
  withProps(({ data }) => ({
    currentProject: data.projectById,
    currentStep: data.stepById,
  })),
  withProps(({ currentProject }) => ({
    steps: currentProject.stepsByProjectId.nodes,
  })),
  branch(({ currentProject, currentStep }) => currentProject == null || currentStep == null, renderNothing),
  withProps(({ currentStep }) => ({
    gboxId: currentStep.gboxByGbox.id,
    meta: currentStep.gboxByGbox.meta,
    frontend: currentStep.gboxByGbox.frontend,
  })),
  connect<any, any, { currentStep: { id: string } }, any>(
    (state) => ({
      componentStates: state.app.componentStates,
    }),
    {
      componentSetState: actionCreators.componentSetState,
    },
  ),
  withProps<any, any>(
    ({
      submitStep,
      uploadFile,
      globalError,
      currentProjectId,
      currentStep,
      rfrProjectStep,
      componentStates,
      componentSetState,
      resetStep,
      gboxId,
      isReport,
      queryBackend,
      isClient,
      steps,
      rfrProjectReport,
      clearErrors,
      interceptStep,
      nextStep,
    }) => ({
      gFuncs: {
        gNextStep: () => {
          nextStep({ reportAtTheEnd: true });
        },
        gInterceptStep: () => {
          interceptStep({ stepId: currentStep.id });
        },
        gSubmitStep: (payload) => {
          submitStep({ ...payload, id: currentStep.id });
        },
        gUpload: (payload) => {
          uploadFile(payload);
        },
        gArgs: currentStep.args,
        gUploadedFiles: currentStep.uploadedFilesByStepId.nodes,
        gImports: currentStep.importsByStepId.nodes.map((e) => ({
          id: e.id,
          injectInto: e.injectInto,
          extractFrom: e.exportByExportId.extractFrom,
          stepRank: e.exportByExportId.stepByStepId.rank,
          gboxTitle: e.exportByExportId.stepByStepId.gboxByGbox.meta.title,
        })),
        gReportGlobalError: globalError,
        gIsReport: isReport,
        gAvailImps: steps
          .filter((x) => x.rank < currentStep.rank)
          .slice()
          .sort((a, b) => a.rank - b.rank)
          .map((x) =>
            x.exportsByStepId.nodes.map((y) => ({
              stepRank: x.rank,
              stepTitle: x.gboxByGbox.meta.title,
              ...y,
            })),
          )
          .reduce((x, y) => x.concat(y || []), []),
        gStatus: currentStep.status,
        gGetState: (path = []) => componentStates.getIn([currentStep.id, ...path]),
        gSetState: (newState, path = []) => componentSetState({ path: [currentStep.id, ...path], newState }),
        gReset: () => resetStep({ stepId: currentStep.id }),
        gResults: currentStep.results,
        gErrors: currentStep.errors,
        gClearErrors: () => clearErrors({ stepId: currentStep.id }),
        gIsClient: isClient,
        /**
         * queries: [
         *   {
         *     id: "c1819f8c9n243",
         *     injectInto: "dim1Data",
         *   },
         *   {
         *     id: "f1hb915t19fr1",
         *     injectInto: "dim2Data",
         *   },
         *   ...
         * ]
         */
        gQueryBackend: (payload) =>
          new Promise((resolve) => {
            queryBackend({
              ...payload,
              cb: (x) => resolve(x),
              audience: gboxId,
            });
          }),
        gGetMultipleExports: (payload) =>
          new Promise((resolve) => {
            queryBackend({
              ...payload,
              cb: (x) => resolve(x),
              endpoint: 'getMultipleExports',
              audience: '__granatum',
            });
          }),
      },
    }),
  ),
);

export default enhance(StepDisplayer);
