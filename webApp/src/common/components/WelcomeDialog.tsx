import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  Typography as T,
  withStyles,
  withWidth,
} from '@material-ui/core';
import React from 'react';
import ReactMarkdown from 'react-markdown';
import { connect } from 'react-redux';
import { compose } from 'recompose';
import actionCreators from '../redux/actionCreators';

const styles = {
  dialogPaper: {
    maxWidth: '100%',
    width: 800,
    // height: '100%',
  },
};

// language=Markdown
const messageMd = `
This is a graphical **single-cell RNA sequencing** (scRNA-Seq) analysis platform for genomics scientists. It has
all the functionality you need for processing, analyzing, and visualizing your scRNA-Seq data.

GranatumX allows for extremely flexible pipeline building. You can build your own analysis pipeline by adding
modules called **G-boxes**. These G-boxes are like lego pieces that you can put together to form a complete study.

GranatumX organizes your project so that you can work on multiple projects, and generates a detailed report with
all the parameters used in each analysis step as well as publication-quality figures and tables
(with legends and captions). It will help you jump start your manuscript writing.

You can start by using our **recommended pipeline recipe** which covers the most common analysis steps appeared in
previous high-impact scRNA-Seq studies. This pipeline will graphically guide you through the analysis of scRNA-seq
data, starting from expression and metadata tables. It has a comprehensive set of modules for

  - preprocessing (normalization, filtering, etc.)
  - visualization (PCA, t-SNE, etc)
  - clustering
  - marker genes identification
  - enrichment analysis
  - cell pseudo-time (estimated trajectory) construction

You can also add steps that perform the following analysis into your pipeline

  - Imputation
  - Interactive outlier removal

Have fun and make discoveries!
`;

const WelcomeDialog: React.FunctionComponent<any> = ({
  currentProjectId,
  updateAddStepDialog,
  updateWelcomeDialog,
  welcomeDialog,
  width,
  classes,
}) => (
  <Dialog
    onClose={() => {
      updateWelcomeDialog({ open: { $set: false } });
    }}
    open={welcomeDialog.open}
    classes={{ paper: classes.dialogPaper }}
    fullScreen={width === 'xs' || width === 'sm'}
  >
    <DialogContent>
      <img
        style={{ margin: '16px auto 32px auto', display: 'block', width: 300, userSelect: 'none' }}
        alt="granatumx_logo_with_text.svg"
        src="/granatumx_logo_with_text.svg"
      />
      <T variant="h6">Welcome to GranatumX!</T>
      <T component={'div' as any}>
        <ReactMarkdown source={messageMd} />
      </T>
    </DialogContent>
    <DialogActions>
      <Button
        onClick={() => {
          updateWelcomeDialog({ open: { $set: false } });
        }}
      >
        Start building your own pipeline
      </Button>
      <Button
        variant="contained"
        color="primary"
        onClick={() => {
          updateWelcomeDialog({ open: { $set: false } });
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
        Use a recommended recipe
      </Button>
    </DialogActions>
  </Dialog>
);

const enhance = compose(
  connect(
    (state: IReduxState) => ({
      currentProjectId: state.app.currentProjectId,
      welcomeDialog: state.app.dialog.welcome,
    }),
    {
      updateAddStepDialog: actionCreators.updateAddStepDialog,
      updateWelcomeDialog: actionCreators.updateWelcomeDialog,
    },
  ),
  withWidth(),
  withStyles(styles),
);

export default enhance(WelcomeDialog);
