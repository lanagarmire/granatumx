import { Button, Fade, Grow, Paper, Typography as T, withStyles } from '@material-ui/core';
import filesize from 'filesize';
import update from 'immutability-helper';
import React from 'react';
import { branch, compose, lifecycle, renderNothing, withHandlers, withProps } from 'recompose';

import { StyleRules } from '@material-ui/core/styles';
import { defaultTo } from '../utils';
import ArgInput from './ArgInput';
import DefaultStepLayout from './DefaultStepLayout';
import LoadingScreen from './LoadingScreen';
import RenderedResult from './RenderedResult';

const styles: StyleRules<string> = {
  title: {
    marginBottom: '8px',
  },
  description: {
    marginBottom: '16px',
  },
  sectionTitle: {
    fontSize: '1.5em',
  },
  paper: {
    marginTop: '16px',
    padding: '16px',
  },
  form: {
    marginBottom: '32px',
    // '& > *:first-child': {
    //   marginTop: '0',
    // },
  },
  page: {
    margin: '16px 0 32px 0',
  },
  buttonRow: {
    margin: 0,
    marginLeft: 'auto',
    padding: 4,
    display: 'flex',
    flexDirection: 'row',
    backgroundColor: 'rgba(250, 250, 250, 0.7)',
    justifyContent: 'flex-end',
  },
  button: {
    margin: '4px',
    position: 'relative',
  },
};

const DonePage = ({
  meta,
  ggUploadedFiles,
  ggArgs,
  ggImports,
  gResults,
  gIsReport,
  gReset,
  classes,
  gNextStep,
  gSetState,
}) => (
  <div className={classes.page}>
    <T variant="h2" gutterBottom>
      {meta.title}
    </T>
    {ggUploadedFiles.map((x) => (
      <T key={x.injectInto}>
        {x.label}:{' '}
        <strong>
          {x.meta.name} ({filesize(x.meta.size)})
        </strong>
      </T>
    ))}
    {ggArgs.map((x) => (
      <T key={x.injectInto}>
        {x.label}: <strong>{JSON.stringify(x.value)}</strong>
      </T>
    ))}
    {ggImports.map((x) => (
      <T key={x.injectInto}>
        {x.label}:{' '}
        <strong>
          {x.extractFrom}{' '}
          <T component={'span' as any} variant="caption" style={{ display: 'unset' }}>
            (from step {x.stepRank + 1}: {x.gboxTitle})
          </T>
        </strong>
      </T>
    ))}
    <div className={classes.form}>
      {gResults ? (
        gResults.map((el, i) => (
          <Paper elevation={1} classes={{ root: classes.paper }} key={i}>
            <RenderedResult {...el} />
          </Paper>
        ))
      ) : (
        <T>There are no results for this step.</T>
      )}
    </div>
    {!gIsReport && (
      <div className={classes.buttonRow}>
        <Button
          color="secondary"
          onClick={() => {
            gSetState(false, ['disableSubmit']);
            gReset();
          }}
          className={classes.button}
        >
          Reset
        </Button>
        <Button variant="contained" color="primary" onClick={gNextStep} className={classes.button}>
          Next step
        </Button>
      </div>
    )}
  </div>
);

const FrontEndFromSpecs = ({
  ggArgs,
  classes,
  frontend,
  gAvailImps,
  gClearErrors,
  gErrors,
  gReportGlobalError,
  gIsReport,
  gNextStep,
  gReset,
  gResults,
  ggUploadedFiles,
  gGetState,
  gStatus,
  gSetState,
  ggImports,
  meta,
  submitStep,
  gInterceptStep,
}) =>
  gStatus === 'IDLE' ? (
    <DefaultStepLayout
      {...{
        gErrors,
        gClearErrors,
        disableSubmit: gGetState(['disableSubmit']),
        submitStep,
        meta,
      }}
    >
      {[
        { type: 'imports', title: 'Imports', idField: 'injectInto' },
        { type: 'uploadedFiles', title: 'Uploaded Files', idField: 'injectInto' },
        { type: 'args', title: 'Arguments', idField: 'injectInto' },
      ].map(
        (section) =>
          frontend[section.type] &&
          frontend[section.type][0] && (
            <div key={section.type}>
              <T variant="h4" className={classes.sectionTitle}>
                {section.title}
              </T>
              <div className={classes.form}>
                {frontend[section.type].map((comp) => (
                  <Paper elevation={1} classes={{ root: classes.paper }} key={comp[section.idField]}>
                    <ArgInput
                      section={section.type}
                      getState={(path = []) => gGetState([section.type, comp.injectInto || comp.extractFrom, ...path])}
                      setState={(x, path = []) =>
                        gSetState(x, [section.type, comp.injectInto || comp.extractFrom, ...path])
                      }
                      defaultState={comp.default}
                      {...{ comp, gAvailImps }}
                    />
                  </Paper>
                ))}
              </div>
            </div>
          ),
      )}
    </DefaultStepLayout>
  ) : gStatus === 'INITIATED' ? (
    <LoadingScreen label="Initiated ..." />
  ) : gStatus === 'INTERCEPTION_REQUESTED' ? (
    <LoadingScreen label="Stopping ..." />
  ) : gStatus === 'RUNNING' ? (
    <LoadingScreen
      label="Running ..."
      actions={[
        {
          label: 'Cancel',
          action: () => {
            gInterceptStep();
          },
        },
      ]}
      message={gErrors}
    />
  ) : gStatus === 'DONE' ? (
    <DonePage
      {...{
        meta,
        ggArgs,
        ggImports,
        ggUploadedFiles,
        gResults,
        gIsReport,
        gReset,
        classes,
        gNextStep,
        gSetState,
      }}
    />
  ) : (
    <T>Unknown Status</T>
  );

const nullToEmptyArray = (x) => (x == null ? [] : x);

const enhance = compose<any, any>(
  withProps((props) =>
    update(props, {
      frontend: {
        imports: defaultTo([]),
        uploadedFiles: defaultTo([]),
        exports: defaultTo([]),
        args: defaultTo([]),
        results: defaultTo([]),
      },
      gResults: defaultTo([]),
    }),
  ),
  branch(({ gGetState }) => gGetState == null, renderNothing),
  withProps(({ frontend, gArgs, gImports, gUploadedFiles }) => ({
    ggUploadedFiles: gUploadedFiles
      ? gUploadedFiles.map((x) => ({ ...x, ...frontend.uploadedFiles.find((y) => y.injectInto === x.injectInto) }))
      : [],
    ggArgs: gArgs ? gArgs.map((x) => ({ ...x, ...frontend.args.find((y) => y.injectInto === x.injectInto) })) : [],
    ggImports: gImports
      ? gImports.map((x) => ({ ...x, ...frontend.imports.find((y) => y.injectInto === x.injectInto) }))
      : [],
  })),
  withHandlers({
    submitStep: ({
      gReportGlobalError,
      gboxId,
      frontend,
      gQueryBackend,
      gGetState,
      gSetState,
      gSubmitStep,
      gUpload,
    }) => async () => {
      if (
        frontend.imports.some((x) => !x.optional && gGetState(['imports', x.injectInto]) == null) ||
        frontend.uploadedFiles.some((x) => !x.optional && gGetState(['uploadedFiles', x.injectInto, 'file']) == null) ||
        frontend.args.some((x) => !x.optional && gGetState('args', x.injectInto) == null)
      ) {
        gReportGlobalError({
          title: 'Cannot submit step',
          content: 'Some required fields are empty!\n\nAll fields are required unless marked "optional".',
        });
        return;
      }

      gSetState(true, ['disableSubmit']);

      const newDynamicItems = await (async () => {
        if (frontend.dynamic == null) {
          return {};
        }

        const { importInjectIntos, argInjectIntos } = frontend.dynamic;

        return gQueryBackend({
          endpoint: 'dynamic',
          imports: importInjectIntos.map((x) => ({
            exportId: gGetState('imports', x),
            injectInto: x,
          })),
          args: argInjectIntos.map((x) => ({
            value: gGetState(['args', x]),
            injectInto: x,
          })),
        });
      })();

      newDynamicItems.results = defaultTo([])(newDynamicItems.results);
      newDynamicItems.exports = defaultTo([])(newDynamicItems.exports);

      const uploadedFiles = await Promise.all(
        frontend.uploadedFiles.map(
          (x) =>
            new Promise((resolve, reject) => {
              if (gGetState(['uploadedFiles', x.injectInto, 'file']) == null) {
                resolve({ ...x, fileId: null, meta: null });
              } else {
                gUpload({
                  file: gGetState(['uploadedFiles', x.injectInto, 'file']),
                  onProgress: (e) => {
                    if (e.lengthComputable) {
                      gSetState(e.loaded, ['uploadedFiles', x.injectInto, 'loaded']);
                      gSetState(e.total, ['uploadedFiles', x.injectInto, 'total']);
                    }
                  },
                  onLoad: () => null,
                  onError: (e) => {
                    reject(e);
                  },
                  onResponse: (response) => {
                    const { fileId } = JSON.parse(response);
                    gSetState(true, ['uploadedFiles', x.injectInto, 'done']);
                    const fileObj = gGetState(['uploadedFiles', x.injectInto, 'file']);
                    const { name, lastModified, lastModifiedData, size, type } = fileObj;
                    resolve({ ...x, fileId, meta: { name, lastModified, lastModifiedData, size, type } });
                  },
                });
              }
            }),
        ),
      );

      gSubmitStep({
        gbox: gboxId,
        uploadedFiles,
        args: frontend.args.map((x) => ({
          injectInto: x.injectInto,
          value: gGetState(['args', x.injectInto]),
        })),
        imports: frontend.imports.map((x) => ({
          injectInto: x.injectInto,
          exportId: gGetState(['imports', x.injectInto]),
        })),
        exports: [
          ...frontend.exports.map((x) => ({
            extractFrom: x.extractFrom,
            kind: x.kind,
            meta: {},
          })),
          ...newDynamicItems.exports,
        ],
        results: [...frontend.results, ...newDynamicItems.results],
      });

      gSetState(false, ['disableSubmit']);
    },
  }),
  withStyles(styles),
);

export default enhance(FrontEndFromSpecs);
