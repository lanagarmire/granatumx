"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.default = {
    clearRunningTasks: !!process.env.CLEAR_RUNNING_TASKS,
    uploadPathOnHost: '/var/granatum/uploaded_files/',
    dataPathOnHost: '/var/granatum/data/',
    stepPathBaseOnHost: '/var/granatum/steps',
    argsFileOnSWD: 'args.json',
    exportsDirOnSWD: 'exports',
    importsDirOnSWD: 'imports',
    debugDirOnSWD: 'debug',
    uploadedFilesDirOnSWD: 'uploaded_files',
    resultsFileOnSWD: 'results.json',
    exportsAnnoFileOnSWD: 'exports_anno.json',
    // TODO: `errors.json` for reporting errors, and `progress.json` for reporting progress
    swdOnDocker: '/data',
    __DEV__: process.env.NODE_ENV !== 'production',
    GBOX_TIME_LIMIT: process.env.GBOX_TIME_LIMIT != null ? +process.env.GBOX_TIME_LIMIT : 14 * 24 * 60 * 60 * 1000,
};
