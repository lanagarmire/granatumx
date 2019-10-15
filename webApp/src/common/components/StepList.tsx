import { List, Typography as T, withStyles } from '@material-ui/core';
import { common, red } from '@material-ui/core/colors';
import gql from 'graphql-tag';
import _ from 'lodash';
import React from 'react';
import { graphql } from 'react-apollo';
import { DragDropContext, Draggable, Droppable } from 'react-beautiful-dnd';
import { connect } from 'react-redux';
import { branch, compose, renderNothing, withHandlers, withProps } from 'recompose';

import actionCreators from '../redux/actionCreators';
import { guardEmpty, withStateUpdater } from '../utils';
import Scrollbars from './Scrollbars';
import StepListItem from './StepListItem';

const styles = {
  badge: {
    zIndex: 60,
  },
};

const StepList: React.FunctionComponent<any> = ({
  updateAddStepDialog,
  currentStepId,
  currentProjectId,
  classes,
  onScrollToStep,
  steps,
  dragEndHandler,
  rfrProjectStep,
  keyMap,
  prevStep,
  nextStep,
  addStep,
}) => (
  <div style={{ overflowY: 'auto', flex: 1 }}>
    <DragDropContext onDragEnd={dragEndHandler}>
      <Scrollbars>
        {steps.length === 0 ? (
          <div style={{ position: 'relative', height: '100%', width: '100%' }}>
            <div style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%,-50%)' }}>
              <T>
                You don&apos;t seem to have any steps in this project. Try{' '}
                <a
                  role="button"
                  href="#"
                  onClick={() => {
                    updateAddStepDialog({
                      $set: {
                        open: true,
                        showSteps: false,
                        showRecipes: true,
                        projectId: currentProjectId,
                        atStepRank: 0,
                      },
                    });
                  }}
                >
                  use a recipe
                </a>{' '}
                or{' '}
                <a
                  role="button"
                  href="#"
                  onClick={() => {
                    addStep({
                      gboxId: 'UploadFiles',
                      projectId: currentProjectId,
                      atStepRank: 0,
                      jumpToNewStep: true,
                    });
                  }}
                >
                  add an upload step
                </a>
                {/*<a
                    role="button"
                    href="#"
                    onClick={() => {
                      updateAddStepDialog({
                        $set: {
                          open: true,
                          showSteps: true,
                          showRecipes: false,
                          projectId: currentProjectId,
                          atStepRank: 0,
                        },
                      });
                    }}
                  >
                    add a step
                  </a>*/}
              </T>
            </div>
          </div>
        ) : (
          <Droppable droppableId="droppable">
            {({ droppableProps, innerRef: droppableRef, placeholder: droppablePlaceholder }, { isDraggingOver }) => (
              <div ref={droppableRef} {...droppableProps}>
                <List style={{ overflow: 'hidden', paddingBottom: 36 }}>
                  <div style={{ padding: '6px 16px 16px 16px' }}>
                    <T variant="h5">Steps</T>
                  </div>
                  {steps.map((step, index) => (
                    <Draggable key={step.id} draggableId={step.id} index={index}>
                      {(
                        { draggableProps, innerRef: draggableRef, placeholder: draggablePlaceholder, dragHandleProps },
                        { isDragging },
                      ) => (
                        <div>
                          <div
                            ref={draggableRef}
                            {...draggableProps}
                            className={isDragging ? classes.dragging : ''}
                            style={{ display: 'flex', ...draggableProps.style }}
                          >
                            <StepListItem
                              {...{
                                step,
                                dragHandleProps,
                                rfrProjectStep,
                                currentProjectId,
                                currentStepId,
                                onScrollToStep,
                              }}
                            />
                          </div>
                          {isDraggingOver ? (
                            draggablePlaceholder
                          ) : isDragging ? (
                            <div
                              style={{
                                width: '100%',
                                height: 30,
                                backgroundColor: red[300],
                                position: 'relative',
                              }}
                            >
                              <T
                                variant="caption"
                                style={{
                                  color: common.white,
                                  position: 'absolute',
                                  top: '50%',
                                  left: '50%',
                                  transform: 'translate(-50%, -50%)',
                                }}
                              >
                                Remove step
                              </T>
                            </div>
                          ) : null}
                        </div>
                      )}
                    </Draggable>
                  ))}
                </List>
                {droppablePlaceholder}
              </div>
            )}
          </Droppable>
        )}
      </Scrollbars>
    </DragDropContext>
  </div>
);

const enhance: any = compose(
  connect(
    (state: IReduxState) => ({
      currentProjectId: state.app.currentProjectId,
      currentStepId: state.app.currentStepId,
    }),
    {
      rfrProjectStep: actionCreators.rfrProjectStep,
      rfrProjectReport: actionCreators.rfrProjectReport,
      rfrNoSteps: actionCreators.rfrNoSteps,
      removeStep: actionCreators.removeStep,
      reorderStep: actionCreators.reorderStep,
      updateAddStepDialog: actionCreators.updateAddStepDialog,
      addStep: actionCreators.addStep,
    },
  ),
  branch((props: any) => props.currentProjectId == null, renderNothing),
  graphql(
    gql`
      query StepList($currentProjectId: UUID!) {
        projectById(id: $currentProjectId) {
          id
          rank
          stepsByProjectId {
            nodes {
              id
              rank
              status
              errors
              gboxByGbox {
                id
                meta
              }
            }
          }
        }
      }
    `,
  ),
  branch(({ data }) => data.loading, renderNothing),
  guardEmpty('projectById'),
  branch(({ data }) => data.projectById == null, renderNothing),
  withProps(({ data }) => ({
    steps: _(data.projectById.stepsByProjectId.nodes)
      .orderBy('rank')
      .value(),
  })),
  branch(({ steps }) => steps == null, renderNothing),
  withStateUpdater('needsToScrollToStep', true),
  withHandlers({
    onScrollToStep: ({ setNeedsToScrollToStep, needsToScrollToStep }) => ({ stepElement }) => {
      if (needsToScrollToStep) {
        setNeedsToScrollToStep(false);
        setTimeout(() => {
          stepElement.scrollIntoView({ behavior: 'smooth' });
        }, 200);
      }
    },
    dragEndHandler: ({ currentProjectId, currentStepId, rfrNoSteps, reorderStep, removeStep }) => (e) => {
      if (e.destination === null) {
        removeStep({
          stepId: e.draggableId,
          doRedirect: {
            currentStepId,
            currentProjectId,
          },
        });
        return;
      }

      reorderStep({
        stepId: e.draggableId,
        toRank: e.destination.index,
        cb: () => {
          if (currentStepId == null) {
            rfrNoSteps({ projectId: currentProjectId });
          }
        },
      });
    },
  }),
  withProps(({ steps, currentStepId }) => ({
    currentStep: steps.find((s) => s.id === currentStepId),
  })),
  withStyles(styles),
);

export default enhance(StepList);
