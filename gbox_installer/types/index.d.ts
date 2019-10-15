interface IPackageSpec {
  id: string;
  meta?: {
    maintainer?: {
      name: string;
      email?: string;
    };
  };
  gboxes: Array<string | IGboxSpec>;
  buildCommand?: string;
}

interface IRecipeSpec {
  id: string;
  meta: {
    title: string;
    subtitle: string;
    featured?: boolean;
    maintainer?: {
      name: string;
      email?: string;
    };
  };
  steps: Array<{
    gbox: string;
    initialState: any;
  }>;
}

interface IGboxSpec {
  id: string;
  meta: {
    title: string;
    subtitle: string;
    description?: string;
    maintainer?: {
      name: string;
      email?: string;
    };
    citations?: Array<{
      bibtex: string;
    }>;
  };
  endpoints?: object;
  // endpoints?: {
  //   [EndpointIdentifier in string]:
  //     | string
  //     | {
  //         type: 'docker';
  //         image: string;
  //         cmd: string;
  //       }
  // };
  frontend:
    | string
    | {
        uploadedFiles?: Array<{
          injectInto: string;
          optional?: boolean;
          label: string;
          description?: string;
        }>;
        args?: Array<{
          type: 'number' | 'text' | 'seed' | 'select' | 'checkbox' | 'samplePicker';
          injectInto: string;
          default: any;
          label: string;
          description?: string;
          choices?: Array<{
            value: any;
            label?: string;
            description?: string;
          }>;
          min?: number;
          max?: number;
          step?: number;
        }>;
        imports?: Array<{
          kind: 'assay' | 'geneMeta' | 'sampleMeta' | 'sampleCoords' | 'geneCoords';
          injectInto: string;
          label: string;
          description?: string;
        }>;
        exports?: Array<{
          kind: 'assay' | 'geneMeta' | 'sampleMeta' | 'sampleCoords' | 'geneCoords';
          extractFrom: string;
          meta?: any;
        }>;
      };
}
