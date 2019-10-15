import { GraphQLBoolean, GraphQLInt, GraphQLNonNull, GraphQLObjectType, GraphQLString } from 'graphql';
import { Plugin } from 'graphile-build';

export default (knex) => {
  const confirmAccessToProject = async (projectId, userId) => {
    const projectMaybe = (await knex('project')
      .select('owner_id')
      .where({ id: projectId }))[0];

    if (projectMaybe == null) {
      throw new Error('No such project.');
    }

    const ownerId = projectMaybe.owner_id;

    if (ownerId !== userId) {
      throw new Error('User id incorrect. Permission denied.');
    }
  };

  const plugins: Plugin[] = [
    (builder, { pgExtendedTypes }) => {
      builder.hook('GraphQLObjectType:fields', (
        fields, // Input object - the fields for this GraphQLObjectType
        { getTypeByName, newWithHooks }, // Build object - handy utils
        { scope: { isRootQuery, isRootMutation } }, // Context object - used for filtering
      ) => {
        const JSONType = getTypeByName('JSON');
        const UUID = getTypeByName('UUID');

        if (isRootQuery) {
          // const { whoami, nodeId } = fields;
          // return { whoami, nodeId };
          return fields;
        }

        if (isRootMutation) {
          return {
            ...fields,
            // saveStep: {
            //   type: newWithHooks(
            //     GraphQLObjectType,
            //     {
            //       name: 'SaveStepType',
            //       fields: {
            //         success: { type: GraphQLBoolean },
            //       },
            //     },
            //     { isMutationPayload: true },
            //   ),
            //   args: {
            //     input: newWithHooks(
            //       GraphQLObjectType,
            //       {
            //         name: 'SaveStepInput',
            //         fields: {
            //           stepId: { type: new GraphQLNonNull(UUID) },
            //           args: { type: JSONType },
            //           imports: { type: JSONType },
            //           exports: { type: JSONType },
            //           results: { type: JSONType },
            //         },
            //       },
            //       { isMutationInput: true },
            //     ),
            //   },
            //   resolve: async (data, { stepId, args, imports, exports, results }, { user }) => {
            //     console.log('args =', args);
            //     console.log('results =', results);
            //     console.log('imports =', imports);
            //     const step = (await knex('step')
            //       .select('status')
            //       .where({ id: stepId }))[0];
            //
            //     if (step == null) {
            //       throw new Error(`Step with id ${stepId} not found.`);
            //     }
            //
            //     if (step.status !== 'idle') {
            //       throw new Error(`Step is not idle.`);
            //     }
            //
            //     await knex.transaction(async trx => {
            //       await trx('step')
            //         .update({ args: JSON.stringify(args), results: JSON.stringify(results) })
            //         .where({ id: stepId });
            //
            //       await trx('import').insert(
            //         imports
            //           ? imports.map(({ exportId, injectInto }) => ({
            //               step_id: stepId,
            //               export_id: exportId,
            //               inject_into: injectInto,
            //             }))
            //           : [],
            //       );
            //
            //       await trx('export').insert(
            //         exports
            //           ? exports.map(({ kind, name, meta, extractFrom, data }) => ({
            //               step_id: stepId,
            //               kind,
            //               name,
            //               meta,
            //               extract_from: extractFrom,
            //               data,
            //             }))
            //           : [],
            //       );
            //     });
            //
            //     return { success: true };
            //   },
            // },
            modifyProject: {
              type: newWithHooks(
                GraphQLObjectType,
                {
                  name: 'ModifyProjectReturnType',
                  fields: {
                    success: { type: GraphQLBoolean },
                  },
                },
                { isMutationPayload: true },
              ),
              args: {
                id: { type: new GraphQLNonNull(UUID) },
                rank: { type: GraphQLInt },
                name: { type: GraphQLString },
                description: { type: GraphQLString },
              },
              resolve: async (data, { id, rank, name, description }, { user }) => {
                await confirmAccessToProject(id, user.user_id);

                if (rank == null && name == null && description == null) {
                  throw new Error('At least one of rank, name, and description needs to be specified.');
                }

                await knex('project')
                  .where({ id })
                  .update({ rank, name, description })
                  .debug();

                return { success: true };
              },
            },
            // deleteProject: {
            //   type: newWithHooks(
            //     GraphQLObjectType,
            //     {
            //       name: 'DeleteProjectReturnType',
            //       fields: {
            //         success: { type: GraphQLBoolean },
            //       },
            //     },
            //     { isMutationPayload: true },
            //   ),
            //   args: {
            //     id: { type: new GraphQLNonNull(UUID) },
            //   },
            //   resolve: async (data, { id }, { user }) => {
            //     await confirmAccessToProject(id, user.user_id);
            //
            //     await knex('project')
            //       .where({ id })
            //       .del();
            //
            //     return { success: true };
            //   },
            // },
            // createProject: {
            //   type: newWithHooks(
            //     GraphQLObjectType,
            //     {
            //       name: 'CreateProjectReturnType',
            //       fields: {
            //         id: { type: UUID },
            //         rank: { type: GraphQLInt },
            //       },
            //     },
            //     { isMutationPayload: true },
            //   ),
            //   args: {
            //     name: { type: GraphQLString },
            //     description: { type: GraphQLString },
            //   },
            //   resolve: async (data, { name, description }, { user }) => {
            //     let newProject;
            //
            //     await knex.transaction(async (trx) => {
            //       const count = (await trx('project')
            //         .count()
            //         .where({ owner_id: user.user_id }))[0].count;
            //
            //       newProject = (await trx('project')
            //         .insert({ owner_id: user.user_id, rank: count, name, description })
            //         .returning('*'))[0];
            //
            //       const defaultGboxes = [
            //         'UploadFiles',
            //         'Imputation',
            //         'DimReduction',
            //         'InteractiveOutlierRemoval',
            //         'DimReduction',
            //         'Normalization',
            //         'GeneFiltering',
            //         'SampleClustering',
            //         'DifferentialExpression',
            //         'MonoclePseudotime',
            //       ];
            //
            //       const steps = defaultGboxes.map((x, i) => ({ project_id: newProject.id, rank: i, gbox: x }));
            //
            //       await trx('step').insert(steps);
            //     });
            //
            //     return newProject;
            //   },
            // },
          };
        }

        return fields;
      });
    },
  ];

  return plugins;
};
