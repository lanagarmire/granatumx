import {
  Button,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  ExpansionPanel,
  ExpansionPanelDetails,
  ExpansionPanelSummary,
  List,
  Menu,
  MenuItem,
  Tab,
  Tabs,
  Tooltip,
  Typography as T,
  withStyles,
} from '@material-ui/core';
import { blue, cyan, green, orange } from '@material-ui/core/colors';
import { ExpandMore } from '@material-ui/icons';
import gql from 'graphql-tag';
import update from 'immutability-helper';
import _ from 'lodash';
import React from 'react';
import { graphql } from 'react-apollo';
import JSONTree from 'react-json-tree';
import { connect } from 'react-redux';
import { branch, renderNothing, withHandlers, withProps, withState } from 'recompose';
import { compose } from 'redux';

import { QUERY_BACKEND } from '../constants';
import actionCreators from '../redux/actionCreators';
import { downloadURI, guardEmpty } from '../utils';
import LoadingScreen from './LoadingScreen';
import RenderedResult from './RenderedResult';

const styles: any = (theme) => ({
  chip: {
    margin: theme.spacing.unit * 0.75,
    cursor: 'pointer',
  },
  chipLabel: {
    display: 'block',
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    maxWidth: 220,
  },
  chip_INITIATED: {
    opacity: 0.3,
  },
  chip_RUNNING: {
    opacity: 0.3,
  },
  chip_sampleCoords: {
    backgroundColor: cyan[100],
  },
  chip_assay: {
    backgroundColor: blue[100],
  },
  chip_geneMeta: {
    backgroundColor: orange[100],
  },
  chip_sampleMeta: {
    backgroundColor: green[100],
  },
  group: {
    margin: -theme.spacing.unit * 0.75,
    display: 'flex',
    justifyContent: 'left',
    flexWrap: 'wrap',
    paddingTop: 0,
  },
  tooltipTitle: {
    color: '#fff',
  },
  previewDialog: {
    width: 800,
    height: 800,
    maxWidth: '90vw',
    maxHeight: '90vh',
  },
  metaTable: {
    width: 'auto',
    margin: 'auto',
    marginTop: 16,
  },
  metaTableRow: {
    height: 32,
  },
  assayHeatMapCell: {
    width: 10,
    height: 10,
  },
});

const DataManager = ({
  queryBackend,
  previewDialog,
  updatePreviewDialog,
  setQuerying,
  querying,
  exportsGroupedByStep,
  folding,
  updateFolding,
  classes,
  currentStepId,
  downloadGPiece,
  anchor,
  setAnchor,
  globalDialog,
}) => (
  <>
    <List>
      <div style={{ padding: 16 }}>
        <T variant="h5">Project data</T>
      </div>
      {exportsGroupedByStep.length === 0 ? (
        <div style={{ margin: 16 }}>
          <T variant="caption">
            You currently don&apos;t have any data in this project. Once you run a step, its exported data will appear
            here.
          </T>
        </div>
      ) : (
        exportsGroupedByStep.map((x, i) => (
          <ExpansionPanel
            key={x.id}
            expanded={folding[i] || x.id === currentStepId}
            onChange={(_e, b) => updateFolding({ [i]: { $set: b } })}
          >
            <ExpansionPanelSummary expandIcon={x.id === currentStepId ? null : <ExpandMore />}>
              <T>
                Step {x.rank + 1}: {x.gboxByGbox.meta.title}
              </T>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails classes={{ root: classes.group }}>
              {x.exports.map((piece) => (
                <Tooltip key={piece.id} title={piece.kind}>
                  <Chip
                    className={`${classes.chip} ${classes[`chip_${piece.kind}`]} ${classes[`chip_${piece.status}`]}`}
                    classes={{ label: classes.chipLabel }}
                    label={piece.extractFrom}
                    onClick={({ target }) => {
                      if (piece.status === 'DONE') {
                        setAnchor({ target, piece });
                      }
                    }}
                  />
                </Tooltip>
              ))}
            </ExpansionPanelDetails>
          </ExpansionPanel>
        ))
      )}
    </List>
    <Menu
      id="simple-menu"
      anchorEl={anchor.target}
      open={!!anchor.target}
      onClose={() => {
        setAnchor({ target: null, piece: null });
      }}
    >
      {/* <MenuItem
        disabled={querying}
        style={{ position: 'relative' }}
        onClick={() => {
          setQuerying(true);
          queryBackend({
            audience: '__granatum',
            endpoint: 'getDataPiecePreview',
            cb: (result) => {
              updatePreviewDialog({ $merge: { piece: anchor.piece, ...result } });
              setQuerying(false);
            },
            piece: {
              extractFrom: anchor.piece.extractFrom,
              exportId: anchor.piece.id,
            },
          });
          setAnchor({ target: null, piece: null });
          updatePreviewDialog({ open: { $set: true } });
        }}
      >
        Preview
        {querying && (
          <CircularProgress
            style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)' }}
            size={24}
          />
        )}
      </MenuItem> */}
      {__DEV__ && (
        <MenuItem
          onClick={() => {
            globalDialog({
              content: <JSONTree data={anchor.piece} />,
            });
          }}
        >
          Info
        </MenuItem>
      )}
      <MenuItem
        onClick={() => {
          if (anchor.piece != null) {
            downloadGPiece(anchor.piece);
            setAnchor({ target: null, piece: null });
          }
        }}
      >
        Download
      </MenuItem>
    </Menu>
    <Dialog
      classes={{ paper: classes.previewDialog }}
      open={previewDialog.open}
      onClose={() => {
        updatePreviewDialog({ open: { $set: false } });
      }}
    >
      {querying ? (
        <DialogContent>
          <LoadingScreen />
        </DialogContent>
      ) : (
        <>
          {previewDialog.extractFrom && <DialogTitle>{previewDialog.extractFrom}</DialogTitle>}
          <Tabs
            value={previewDialog.show}
            onChange={(e, v) => updatePreviewDialog({ show: { $set: v } })}
            indicatorColor="primary"
            textColor="primary"
            centered
          >
            <Tab label="Rendered" value="rendered" />
            <Tab label="Source" value="source" />
          </Tabs>
          <DialogContent>
            {previewDialog.show === 'rendered' ? (
              previewDialog.rendered ? (
                <RenderedResult {...previewDialog.rendered} />
              ) : (
                <div style={{ display: 'flex' }}>
                  <T
                    variant="caption"
                    component={'div' as any}
                    style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50% ,-50%)' }}
                  >
                    <p>
                      We don&apos;t know how to render data with <strong>kind: {previewDialog.kind}</strong>
                    </p>
                    <p>
                      You might want to take a look at its{' '}
                      <a
                        role="button"
                        href="#"
                        onClick={(e) => {
                          updatePreviewDialog({ show: { $set: 'source' } });
                          e.preventDefault();
                        }}
                      >
                        source
                      </a>
                      .
                    </p>
                  </T>
                </div>
              )
            ) : previewDialog.show === 'source' ? (
              previewDialog.sourcePreview ? (
                <RenderedResult {...previewDialog.sourcePreview} />
              ) : (
                <div style={{ display: 'flex' }}>
                  <T
                    variant="caption"
                    component={'div' as any}
                    style={{ position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50% ,-50%)' }}
                  >
                    <p>Server did not return any source preview.</p>
                  </T>
                </div>
              )
            ) : null}
          </DialogContent>
        </>
      )}

      <DialogActions>
        <Button
          onClick={() => {
            if (previewDialog.piece != null) {
              downloadGPiece(previewDialog.piece);
            }
          }}
        >
          Download
        </Button>
        <Button
          color="primary"
          onClick={() => {
            updatePreviewDialog({ open: { $set: false } });
          }}
        >
          Dismiss
        </Button>
      </DialogActions>
    </Dialog>
  </>
);

const enhance = compose(
  connect(
    (state: IReduxState) => ({
      currentProjectId: state.app.currentProjectId,
      currentStepId: state.app.currentStepId,
    }),
    {
      queryBackend: actionCreators.queryBackend,
      globalDialog: actionCreators.globalDialog,
    },
  ),
  branch((props: any) => props.currentProjectId == null, renderNothing),
  withHandlers({
    downloadGPiece: () => (piece) => {
      downloadURI(`/download-data/${piece.id}`);
      // queryBackend({
      //   audience: '__granatum',
      //   endpoint: 'getDataDownloadLink',
      //   piece: { exportId: piece.id },
      //   cb: ({ fileId, filename }) => {
      //     downloadURI(`/download-file/${fileId}`, filename);
      //   },
      // });
    },
  }),
  graphql(gql`
    query DataManager($currentProjectId: UUID!) {
      projectById(id: $currentProjectId) {
        id
        stepsByProjectId {
          nodes {
            id
            status
            exportsByStepId {
              nodes {
                id
                kind
                meta
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
        }
      }
    }
  `),
  branch(({ data }) => data.loading, renderNothing),
  guardEmpty('projectById'),
  branch(({ data }) => data.projectById == null, renderNothing),
  withProps(({ data }) => ({
    exportsGroupedByStep: _(data.projectById.stepsByProjectId.nodes)
      .flatMap((x) => x.exportsByStepId.nodes.map((y) => ({ ...y, status: x.status })))
      .groupBy((x) => x.stepByStepId.id)
      .map((x: any) => ({ ...x[0].stepByStepId, exports: _.sortBy(x, ['kind', 'extractFrom']) }))
      .sortBy('rank')
      .value(),
  })),
  withState('contextMenu', 'setContextMenu', false),
  withState('anchor', 'setAnchor', { target: null, piece: null }),
  withState('querying', 'setQuerying', false),
  withState('previewDialog', 'setPreviewDialog', {
    open: false,
    piece: null,
    show: 'rendered',
    extractFrom: null,
    data: null,
  }),
  withState('folding', 'setFolding', []),
  withHandlers({
    updateFolding: ({ folding, setFolding }) => (updater) => setFolding(update(folding, updater)),
    updatePreviewDialog: ({ previewDialog, setPreviewDialog }) => (updater) =>
      setPreviewDialog(update(previewDialog, updater)),
  }),
  withStyles(styles),
);

export default enhance(DataManager);
