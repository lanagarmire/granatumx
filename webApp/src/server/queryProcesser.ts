import cp from 'child_process';
import { readJsonSync, writeFileSync } from 'fs-extra';
import _ from 'lodash';
import Papa from 'papaparse';
import { resolve } from 'path';
import uuidv4 from 'uuid/v4';

import { trace } from '../common/utils';
import config from './config';

import getSingletonKnex from './getSingletonKnex';

// TODO: Security issue: the client could potentially access other people's data

let knex;

const getInput = async (args, imports) => {
  const input = {};

  if (args) {
    args.forEach((a) => {
      input[a.injectInto] = a.value;
    });
  }

  if (imports) {
    await Promise.all(
      imports.map(async (d) => {
        input[d.injectInto] = readJsonSync(resolve(config.dataPath, d.exportId));
      }),
    );
  }

  return input;
};

const getDataPiece = async (piece) => {
  if (piece == null) {
    throw new Error(`The "piece" is missing.`);
  }

  return {
    ...(await knex('export')
      .select(['kind', 'name'])
      .where('id', piece.exportId))[0],
    data: readJsonSync(resolve(config.dataPath, piece.exportId)),
  };
};

/**
 * query.piece = { name: 'My Assay', exportId: '8fe8a-9fea..' },
 */
const getDataPiecePreview = async ({ piece }) => {
  const exp = await getDataPiece(piece);

  const rendered = (() => {
    switch (exp.kind) {
      case 'assay':
        return {
          type: 'markdown',
          data: [
            `Number of cells: **${exp.data.sampleIds.length}**`,
            '',
            '```',
            exp.data.sampleIds.slice(0, 10).join(', '),
            '```',
            '',
            `Number of genes: **${exp.data.geneIds.length}**`,
            '',
            '```',
            exp.data.geneIds.slice(0, 10).join(', '),
            '```',
            '',
            'Top-left corner of the matrix:',
            '',
            '```',
            exp.data.matrix
              .slice(0, 10)
              .map((x) => x.slice(0, 10).join(', '))
              .join('\n'),
            '```',
          ].join('\n'),
          // data: {
          //   sampleIds: exp.data.sampleIds,
          //   geneIds: exp.data.geneIds.slice(0, 100),
          //   matrix: exp.data.matrix.slice(0, 100),
          // },
          height: 600,
          width: 750,
        };
      case 'sampleMeta':
        return {
          type: 'compactTable',
          data: {
            cols: ['Sample ID', 'Value'],
            data: Object.entries(exp.data).map(([k, v]) => ({ 'Sample ID': k, Value: v })),
          },
        };
      case 'geneMeta':
        return {
          type: 'compactTable',
          data: {
            cols: ['Gene ID', 'Value'],
            data: Object.entries(exp.data).map(([k, v]) => ({ 'Gene ID': k, Value: v })),
          },
        };
      case 'sampleCoords':
        return {
          type: 'scatterplot',
          height: 600,
          width: 750,
          data: exp.data,
        };
      default:
        return null;
    }
  })();

  const source = JSON.stringify(exp);

  const LENGTH_LIMIT = 10000;

  const sourcePreview = (() => {
    if (source.length > LENGTH_LIMIT + 50) {
      return {
        type: 'codeExcerpt',
        data: {
          head: source.slice(0, Math.floor(LENGTH_LIMIT / 2)),
          tail: source.slice(source.length - Math.floor(LENGTH_LIMIT / 2)),
          numBytesOmitted: source.length - LENGTH_LIMIT,
        },
      };
    }

    return {
      type: 'code',
      data: source,
    };
  })();

  return trace('queryProcesser is returning', {
    ...piece,
    rendered,
    sourcePreview,
  });
};

const formatPiece = (exp) => {
  switch (exp.kind) {
    case 'assay': {
      const table = [
        ['Gene ID', ...exp.data.sampleIds],
        ..._.zipWith(exp.data.geneIds, exp.data.matrix, (x, y: any) => [x, ...y]),
      ];
      return { filename: `${exp.name}.csv`, content: Papa.unparse(table) };
    }
    case 'sampleMeta': {
      return {
        filename: `${exp.name}.csv`,
        content: Papa.unparse(
          _(exp.data)
            .map((v, k) => ({ 'Sample ID': k, [exp.name]: v }))
            .value(),
        ),
      };
    }
    case 'geneMeta': {
      return {
        filename: `${exp.name}.csv`,
        content: Papa.unparse(
          _(exp.data)
            .map((v, k) => ({ 'Gene ID': k, [exp.name]: v }))
            .value(),
        ),
      };
    }
    case 'raw': {
      return {
        filename: exp.name,
        content: exp.data,
      };
    }
    default: {
      return {
        filename: `${exp.name}.json`,
        content: JSON.stringify(exp.data),
      };
    }
  }
};

const getMultipleExports = async ({ exportIds }, user) =>
  Promise.all(
    exportIds.map(async (exportId) => {
      const res = await knex('export')
        .join('step', 'step.id', 'export.step_id')
        .join('project', 'project.id', 'step.project_id')
        .select(['project.owner_id'])
        .where({ 'export.id': exportId });

      const ownerId = res[0].owner_id;

      if (user.user_id !== ownerId) {
        throw new Error(`You do not have access to this export.`);
      }

      return readJsonSync(resolve(config.dataPath, exportId));
    }),
  );

const getDataDownloadLink = async ({ piece }, user) => {
  const exp = await getDataPiece(piece);
  const { filename, content } = formatPiece(exp);
  const fileId = uuidv4();

  writeFileSync(resolve(config.downloadFilePath, fileId), content);

  return { fileId, filename };
};

const processModuleQuery = async (query) => {
  const { audience, endpoint, imports, args } = query;

  const exeCommand = (await knex('gbox')
    .select('endpoints')
    .where('id', audience))[0].endpoints[endpoint].split(/ +/);

  const input = await getInput(args, imports);

  const stdin = JSON.stringify(input);

  if (__DEV__) {
    writeFileSync('/tmp/endpoint_input.json', stdin);
  }

  const result = cp.spawnSync(exeCommand[0], exeCommand.slice(1), { input: stdin, cwd: resolve('../pwd') });
  if (result.error != null) {
    throw new Error(
      __DEV__
        ? `message: ${result.error.message}, stack: ${result.error.stack}`
        : 'An internal error just happened. Please try again.',
    );
  }

  const stderr = result.stderr.toString();
  const stdout = result.stdout.toString();

  if (__DEV__) {
    writeFileSync('/tmp/endpoint_output.json', stdout);
    writeFileSync('/tmp/endpoint_stderr', stderr);
  }

  let output;

  try {
    output = JSON.parse(stdout);
  } catch (e) {
    throw new Error(
      __DEV__ ? `message: ${e.message}, stderr: ${stderr}` : 'An internal error just happened. Please try again.',
    );
  }

  return output;
};

const saveStep = async ({ stepId, args, imports, uploadedFiles, exports, results }, user) => {
  const step = (await knex('step')
    .select('status')
    .where({ id: stepId }))[0];

  if (step == null) {
    throw new Error(`Step with id ${stepId} not found.`);
  }

  if (step.status !== 'idle') {
    throw new Error(`Step is not idle.`);
  }

  await knex.transaction(async (trx) => {
    await trx('step')
      .update({ args: JSON.stringify(args), results: JSON.stringify(results) })
      .where({ id: stepId });

    await trx('import').insert(
      imports
        ? imports.map(({ exportId, injectInto }) => ({
            step_id: stepId,
            export_id: exportId,
            inject_into: injectInto,
          }))
        : [],
    );

    if (exports && exports.length > 0) {
      const exportIds = await trx('export')
        .insert(
          exports
            ? exports.map(({ kind, name, meta, extractFrom }) => ({
                step_id: stepId,
                kind,
                name,
                meta,
                extract_from: extractFrom,
              }))
            : [],
        )
        .returning('id');

      await Promise.all(
        exportIds.map(async (exportId, i) => {
          writeFileSync(resolve(config.dataPath, exportId), JSON.stringify(exports[i].data));
        }),
      );
    }

    await trx('uploaded_file')
      .insert(
        uploadedFiles
          ? uploadedFiles
              .filter((x) => x.fileId != null)
              .map(({ injectInto, fileId, meta }) => ({
                id: fileId,
                step_id: stepId,
                inject_into: injectInto,
                meta,
              }))
          : [],
      )
      .debug();
  });
};

const tryProcessingQuery = async ({ query, user }) => {
  if (query == null) {
    throw new Error('Empty query.');
  }

  const { audience, endpoint } = query;

  if (audience == null) {
    throw new Error('No audience specified missing.');
  }

  if (endpoint == null) {
    throw new Error('Endpoint missing.');
  }

  if (audience === '__granatum') {
    switch (endpoint) {
      case 'getDataPiecePreview':
        return getDataPiecePreview(query);
      case 'getDataDownloadLink':
        return getDataDownloadLink(query, user);
      case 'getMultipleExports':
        return getMultipleExports(query, user);
      case 'saveStep':
        return saveStep(query, user);
      default:
        throw new Error(`Unknown endpoint (${endpoint}) for the audience __granatum.`);
    }
  }

  return processModuleQuery(query);
};

export default async () => {
  knex = await getSingletonKnex();

  return async (req, res) => {
    let response;

    try {
      response = await tryProcessingQuery({ query: req.body, user: req.user });
    } catch (e) {
      res.status(500).json({ error: e.message });
      console.error(e);
      return;
    }

    if (response == null) {
      response = null;
    }

    res.json(response);
  };
};
