import chalk from 'chalk';
import {
  ensureDirSync,
  ensureFileSync,
  ensureLinkSync,
  outputJsonSync,
  pathExistsSync,
  readJsonSync,
  removeSync,
} from 'fs-extra';
import Knex from 'knex';
import _ from 'lodash';
import { join, resolve } from 'path';
import config from './config';
import { getExeCmdArgv } from './getExeCmd';
import runGbox from './runGbox';
import { safeReadFileSync, timeout } from './utils';

export default async (knex: Knex) => {
  // noinspection InfiniteLoopJS
  while (true) {
    // const stepToRun = (await knex('step')
    //   .select(['id', 'gbox', 'args', 'results', 'project_id'])
    //   .where('status', 'initiated')
    //   .orderBy('updated_at', 'asc')
    //   .limit(1))[0];

    /**
     * Here, you need to do an join + aggregated query to get the step_id,
     * and try to skip locked + update it. If it fails, just restart the loop
     *
     * For limiting the number of instances for a single user: create a table "task_runner_instance(id, step_id)"
     * and use some joining magic to determine how many task_runner_instances a certain user has.
     */

    const stepsSatisfyConditions = (await knex.raw(`
    
select
  step.id
from step
  left join import on (import.step_id = step.id)
  left join export on (import.export_id = export.id)
where
  status = 'initiated'
group by step.id
having
  bool_and(export.is_populated) = true or bool_and(export.is_populated) is null
order by step.updated_at asc
limit 1;

        `)).rows;

    if (stepsSatisfyConditions.length === 0) {
      console.log('No task found.');
      await timeout(1000);
      continue;
    }

    console.log('stepsSatisfyConditions =', stepsSatisfyConditions);

    const stepToRunId = stepsSatisfyConditions[0].id;
    console.log('stepToRunId =', stepToRunId);

    let stepToRun: any;

    try {
      stepToRun = await knex.transaction(async (trx) => {
        const res = await trx.raw(
          `
          select
            id,
            gbox,
            args,
            results,
            project_id
          from step
          where
            id = ?
          order by step.updated_at asc
          limit 1
          for update skip locked;
        `,
          [stepToRunId],
        );

        console.log('res =', res);

        if (res.rows[0] == null) {
          throw new Error('The found task is taken by other task runners.');
        }

        await trx.raw(
          `
          update
            step
          set
            status = 'running'
          where
            id = ?;
        `,
          [res.rows[0].id],
        );

        console.log('res.rows[0] =', res.rows[0]);

        return res.rows[0];
      });
    } catch (e) {
      console.log(`Error: ${e.message}`);
      continue;
    }

    console.log('stepToRun =', stepToRun);

    // await knex('step')
    //   .update('status', 'running')
    //   .where('id', stepToRun.id);

    console.log('Found task:', stepToRun.id);

    /**
     *
     * Warehouse   --Hardlink->   host step base  --mount->  docker step base
     *
     * Here is the structure of the "step working directory (SWD)" for a step being run:
     *
     *    <stepId>/
     *        uploaded-files/               <- read only
     *            file_1.csv
     *            file_2.csv
     *        imports/                     <- read only
     *            assayToBeLoaded.json
     *            groupVector.json
     *        exports/
     *            imWrittenByTheGbox.json
     *        args.json                    <- read only
     *        results.json
     *
     * For steps run in Docker, this folder is mounted on /data.
     *
     * In the future if we are going to use Singularity, consider complying to the SCI-F standard.
     *
     *    https://sci-f.github.io/spec-v1
     *
     */

    try {
      const swdOnHost = resolve(config.stepPathBaseOnHost, stepToRun.id);
      removeSync(swdOnHost);

      const taskRelatedFilesOnSWD: Array<{ path: string; opt?: string }> = [];

      /**
       * The following sections prepare the step working directory (SWD):
       *
       *    - args.json
       *    - uploadedFiles/*
       *    - imports/*
       *    - export/
       *    - results.json (empty)
       *
       * For each of the files being prepared, two things are done.
       *
       */

      // --------- args --------

      const argsAsObject = _.fromPairs(
        stepToRun.args.map((x: { injectInto: string; value: any }) => [x.injectInto, x.value]),
      );
      outputJsonSync(resolve(swdOnHost, config.argsFileOnSWD), argsAsObject);
      taskRelatedFilesOnSWD.push({ path: config.argsFileOnSWD, opt: 'ro' });

      // --------- uploadFiles --------

      const currentUploadedFiles: Array<{ id: string; inject_into: string; meta: any }> = await knex('uploaded_file')
        .select(['id', 'inject_into', 'meta'])
        .where('step_id', stepToRun.id);

      const stepUploadedFilesDir = resolve(config.stepPathBaseOnHost, stepToRun.id, 'uploaded-files');
      ensureDirSync(stepUploadedFilesDir);

      taskRelatedFilesOnSWD.push(
        ...(await Promise.all(
          currentUploadedFiles.map(async (upl: { id: string; inject_into: string; meta: { name: string } }) => {
            const fpGbox = resolve(swdOnHost, config.uploadedFilesDirOnSWD, upl.inject_into, upl.meta.name);
            console.log(
              chalk.bold.blue(`linking uploadedFile:  ${resolve(config.uploadPathOnHost, upl.id)}  ->  ${fpGbox}`),
            );
            removeSync(fpGbox);
            ensureLinkSync(resolve(config.uploadPathOnHost, upl.id), fpGbox);
            return { path: join(config.uploadedFilesDirOnSWD, upl.inject_into, upl.meta.name), opt: 'ro' };
          }),
        )),
      );

      // --------- imports --------

      const currentImports: Array<{ id: string; export_id: string; inject_into: string }> = await knex('import')
        .select(['id', 'export_id', 'inject_into'])
        .where('step_id', stepToRun.id);

      taskRelatedFilesOnSWD.push(
        ...(await Promise.all(
          currentImports.map(async (imp: { inject_into: string; export_id: string }) => {
            const fpGbox = resolve(swdOnHost, config.importsDirOnSWD, imp.inject_into);
            console.log(
              chalk.bold.blue(`linking import:  ${resolve(config.dataPathOnHost, imp.export_id)}  ->  ${fpGbox}`),
            );
            removeSync(fpGbox);
            ensureLinkSync(resolve(config.dataPathOnHost, imp.export_id), fpGbox);
            return { path: join(config.importsDirOnSWD, imp.inject_into), opt: 'ro' };
          }),
        )),
      );

      // --------- exports --------

      const exportsDirOnHost = resolve(swdOnHost, config.exportsDirOnSWD);
      ensureDirSync(exportsDirOnHost);
      taskRelatedFilesOnSWD.push({ path: config.exportsDirOnSWD, opt: '' });

      // --------- exports_anno -------

      const exportsAnnoFileOnHost = resolve(swdOnHost, config.exportsAnnoFileOnSWD);
      ensureFileSync(exportsAnnoFileOnHost);
      taskRelatedFilesOnSWD.push({ path: config.exportsAnnoFileOnSWD, opt: '' });

      // --------- results --------

      const resultsFileOnHost = resolve(swdOnHost, config.resultsFileOnSWD);
      ensureFileSync(resultsFileOnHost);
      taskRelatedFilesOnSWD.push({ path: config.resultsFileOnSWD, opt: '' });

      // --------- debug --------

      // the debug folder is mounted as is
      const debugDirOnHost = resolve(swdOnHost, config.debugDirOnSWD);
      ensureDirSync(debugDirOnHost);
      taskRelatedFilesOnSWD.push({ path: config.debugDirOnSWD, opt: '' });

      /**
       * Next, we generate the cmd from the gbox spec
       */

      const exeCmdArgv = await getExeCmdArgv(knex, stepToRun.gbox, taskRelatedFilesOnSWD, swdOnHost);

      /*---------------------------------------------------------*/
      /*---------------------------------------------------------*/
      /*---------------------------------------------------------*/
      const gboxStdoutFileOnHost = resolve(swdOnHost, 'stdout');
      ensureFileSync(gboxStdoutFileOnHost);
      const gboxStderrFileOnHost = resolve(swdOnHost, 'stderr');
      ensureFileSync(gboxStderrFileOnHost);

      const gboxStartTime = process.hrtime();

      await runGbox(
        knex,
        exeCmdArgv,
        gboxStdoutFileOnHost,
        gboxStderrFileOnHost,
        async () =>
          (await knex('step')
            .select(['status'])
            .where({ id: stepToRun.id }))[0].status === 'interception_requested',
        async (obj) => {
          await knex('step')
              .update({'errors': obj})
              .where({ id: stepToRun.id })
        },
      );

      await knex('step')
          .update({'errors': null})
          .where({ id: stepToRun.id })

      const gboxRunningDuration = process.hrtime(gboxStartTime);

      console.log(
        `Gbox has finished running in ${gboxRunningDuration[0]} seconds and ${gboxRunningDuration[1]} nanoseconds`,
      );

      /*---------------------------------------------------------*/
      /*---------------------------------------------------------*/
      /*---------------------------------------------------------*/

      /**
       * The step fails (and thus the stdout and stderr are reported via "step.errors") when one of the following
       * things happen:
       *
       *    1) an statically or dynamically annotated export is not found in form of a file under `exports/`
       *    2) any *.json file in this protocol is not in the JSON format, or not in the structure specified.
       *
       * The fact that stdout and/or stderr is non-empty does *not* result failing of a step.
       *
       */

      /**
       * Now we extract the results from the `results.json`. If the file doesn't exist or it's not a valid JSON,
       * throw an error which will later be repoarted to the front-end (in development mode).
       *
       * New protocol: there's really no reason for static specs of results. Even for the standard modules.
       *
       * To unify the way we process the results. All results are kept untouched before writing to the database.
       * The standard modules need to write their `results.json` in the array (rows) format.
       *
       */
      const results = (() => {
        try {
          return readJsonSync(resultsFileOnHost);
        } catch (e) {
          throw new Error(`Error reading ${resultsFileOnHost}.`);
        }
      })();

      /**
       * There are two types of annotations. The ones appear in the database (in the "exports" table) are the
       * *static* ones. That is, they are there because they are statically specified in the gbox spec file (YAML).
       * The ones appear in the `exports_annotation.json` are the *dynamic* ones -- as they are dynamically
       * generated by the gbox at *run-time*. We need the dynamic annotation because sometimes the number of exports
       * is unknown before running the gbox. For example, the differential expression step might generate a
       * different number of tables, according to the number of distinct groups in the grouping vector.
       *
       * We first list static annotations. For each "export" database row we look for
       * a file in the exports folder with a filename the same as the "extract_from" field in that database row.
       * If the file is found, we store that file into our data warehouse by making the following hardlink:
       *
       *    /var/granatum/steps/<step-id>/exports/aaa  <-->  /var/granatum/data/<data-id>
       *
       * If the file is *not* found, it means that a statically specified export is not generated from the Gbox.
       * This is considered an error and reported to the frontend (via pg "step.errors"). It would also fail the
       * step, and same actions would be taken.
       *
       * Finally, we look for annotation rows in the `exports_annotation.json`. This JSON file uses the same data
       * rows format you get from a run of Knex query. It looks like this:
       *
       *    [
       *      {
       *        extractFrom: 'aaa',
       *        kind: 'assay',
       *        meta: {...}
       *      },
       *      {
       *        extractFrom: 'bbb',
       *        kind: 'sampleMeta',
       *        meta: {...}
       *      },
       *    ]
       *
       * For each of these rows we first look for the corresponding files, and non-existence of a file results in an
       * error just like the static annotation. If the file exists, we proceed to write the row into the database so
       * that we can obtain an ID, which in turn enables us to make the hardlink same as in the static case.
       *
       * If a dynamically annotated export is also found to be statically annotated, an error will be raised.
       *
       */

      const annotatedExportsFromJsonFile = (() => {
        try {
          return readJsonSync(resolve(swdOnHost, config.exportsAnnoFileOnSWD));
        } catch (e) {
          console.log(
            `No dynamic export annotation: ${config.exportsAnnoFileOnSWD} does not exist or is not parsed correctly.`,
          );
          return [];
        }
      })();

      await Promise.all(
        annotatedExportsFromJsonFile.map(async (exp: { extractFrom: string; kind?: string; meta?: any }) => {
          if (exp.extractFrom == null) {
            throw new Error(
              `The field "extractFrom" not found. This field is required in a dynamic export annotation.`,
            );
          }

          await knex('export').insert({
            step_id: stepToRun.id,
            extract_from: exp.extractFrom,
            kind: exp.kind || null,
            meta: exp.meta || null,
          });
        }),
      );

      const annotatedExportsFromDatabase = await knex('export')
        .select(['id', 'extract_from'])
        .where({ step_id: stepToRun.id });

      await Promise.all(
        annotatedExportsFromDatabase.map(async (exp: { id: string; extract_from: string }) => {
          const expFile = resolve(exportsDirOnHost, exp.extract_from);
          if (!pathExistsSync(expFile)) {
            throw new Error(`The annotated export file ${expFile} does not exist.`);
          }
          removeSync(resolve(config.dataPathOnHost, exp.id));
          console.log(chalk.bold.blue(`linking export:  ${expFile}  ->  ${resolve(config.dataPathOnHost, exp.id)}`));
          ensureLinkSync(expFile, resolve(config.dataPathOnHost, exp.id));
          await knex('export')
            .update({ is_populated: true })
            .where({ id: exp.id });
        }),
      );

      /**
       * Finally we writes results and errors into the database
       */

      await knex.transaction(async (trx) => {
        await trx('step')
          .update({
            status: 'done',
            results: JSON.stringify(results),
            state: JSON.stringify({
              executionDuration: { second: gboxRunningDuration[0], nanoSecond: gboxRunningDuration[1] },
            }),
          })
          .where('id', stepToRun.id);
      });
    } catch (e) {
      const swdOnHost = resolve(config.stepPathBaseOnHost, stepToRun.id);
      const gboxStdout = safeReadFileSync(resolve(swdOnHost, 'stdout')).toString();
      const gboxStderr = safeReadFileSync(resolve(swdOnHost, 'stderr')).toString();
      const gboxLogMessage = [
        gboxStdout ? `---stdout---\n${gboxStdout}` : '---stdout is empty---',
        gboxStderr ? `---stderr---\n${gboxStderr}` : '---stderr is empty---',
      ].join('\n');

      await knex.transaction(async (trx) => {
        // TODO: change to reset_step force
        await trx.raw('SELECT reset_step_recursively(?)', [stepToRun.id]);
        // noinspection JSReferencingMutableVariableFromClosure
        await trx('step')
          .update({
            results: null,
            errors: JSON.stringify(
              config.__DEV__
                ? [{ source: 'NodeJS', message: e.stack }, { source: 'Gbox', message: gboxLogMessage }]
                : [{ source: 'Granatum', messsage: 'An internal error just happened. Please try again.' }],
            ),
          })
          .where('id', stepToRun.id);
      });
      console.error(e);
    }
  }
};
