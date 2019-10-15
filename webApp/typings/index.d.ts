declare const __DEV__: boolean;

type IReduxState = any;

interface IStateGlobalDialog {
  type?: 'info' | 'error' | 'warning';
  id: string;
  title: string;
  // if content is a string, it is intepreted as markdown
  content?: string | React.Component;
  actions: Array<{
    label: string;
    description: string;
    action: any;
  }>;
}
