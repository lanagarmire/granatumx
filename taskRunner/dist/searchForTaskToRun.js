"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : new P(function (resolve) { resolve(result.value); }).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = op[0] & 2 ? y["return"] : op[0] ? y["throw"] || ((t = y["return"]) && t.call(y), 0) : y.next) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [op[0] & 2, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
};
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
var _this = this;
Object.defineProperty(exports, "__esModule", { value: true });
var chalk_1 = __importDefault(require("chalk"));
var fs_extra_1 = require("fs-extra");
var lodash_1 = __importDefault(require("lodash"));
var path_1 = require("path");
var config_1 = __importDefault(require("./config"));
var getExeCmd_1 = require("./getExeCmd");
var runGbox_1 = __importDefault(require("./runGbox"));
var utils_1 = require("./utils");
exports.default = (function (knex) { return __awaiter(_this, void 0, void 0, function () {
    var _loop_1;
    var _this = this;
    return __generator(this, function (_a) {
        switch (_a.label) {
            case 0:
                _loop_1 = function () {
                    var stepsSatisfyConditions, stepToRunId, stepToRun, e_1, swdOnHost_1, taskRelatedFilesOnSWD, argsAsObject, currentUploadedFiles, stepUploadedFilesDir, _a, _b, _c, currentImports, _d, _e, _f, exportsDirOnHost_1, exportsAnnoFileOnHost, resultsFileOnHost_1, debugDirOnHost, exeCmdArgv, gboxStdoutFileOnHost, gboxStderrFileOnHost, gboxStartTime, gboxRunningDuration_1, results_1, annotatedExportsFromJsonFile, annotatedExportsFromDatabase, e_2, swdOnHost, gboxStdout, gboxStderr, gboxLogMessage_1;
                    return __generator(this, function (_g) {
                        switch (_g.label) {
                            case 0: return [4 /*yield*/, knex.raw("\n    \nselect\n  step.id\nfrom step\n  left join import on (import.step_id = step.id)\n  left join export on (import.export_id = export.id)\nwhere\n  status = 'initiated'\ngroup by step.id\nhaving\n  bool_and(export.is_populated) = true or bool_and(export.is_populated) is null\norder by step.updated_at asc\nlimit 1;\n\n        ")];
                            case 1:
                                stepsSatisfyConditions = (_g.sent()).rows;
                                if (!(stepsSatisfyConditions.length === 0)) return [3 /*break*/, 3];
                                console.log('No task found.');
                                return [4 /*yield*/, utils_1.timeout(1000)];
                            case 2:
                                _g.sent();
                                return [2 /*return*/, "continue"];
                            case 3:
                                console.log('stepsSatisfyConditions =', stepsSatisfyConditions);
                                stepToRunId = stepsSatisfyConditions[0].id;
                                console.log('stepToRunId =', stepToRunId);
                                _g.label = 4;
                            case 4:
                                _g.trys.push([4, 6, , 7]);
                                return [4 /*yield*/, knex.transaction(function (trx) { return __awaiter(_this, void 0, void 0, function () {
                                        var res;
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0: return [4 /*yield*/, trx.raw("\n          select\n            id,\n            gbox,\n            args,\n            results,\n            project_id\n          from step\n          where\n            id = ?\n          order by step.updated_at asc\n          limit 1\n          for update skip locked;\n        ", [stepToRunId])];
                                                case 1:
                                                    res = _a.sent();
                                                    console.log('res =', res);
                                                    if (res.rows[0] == null) {
                                                        throw new Error('The found task is taken by other task runners.');
                                                    }
                                                    return [4 /*yield*/, trx.raw("\n          update\n            step\n          set\n            status = 'running'\n          where\n            id = ?;\n        ", [res.rows[0].id])];
                                                case 2:
                                                    _a.sent();
                                                    console.log('res.rows[0] =', res.rows[0]);
                                                    return [2 /*return*/, res.rows[0]];
                                            }
                                        });
                                    }); })];
                            case 5:
                                stepToRun = _g.sent();
                                return [3 /*break*/, 7];
                            case 6:
                                e_1 = _g.sent();
                                console.log("Error: " + e_1.message);
                                return [2 /*return*/, "continue"];
                            case 7:
                                console.log('stepToRun =', stepToRun);
                                // await knex('step')
                                //   .update('status', 'running')
                                //   .where('id', stepToRun.id);
                                console.log('Found task:', stepToRun.id);
                                _g.label = 8;
                            case 8:
                                _g.trys.push([8, 20, , 22]);
                                swdOnHost_1 = path_1.resolve(config_1.default.stepPathBaseOnHost, stepToRun.id);
                                fs_extra_1.removeSync(swdOnHost_1);
                                taskRelatedFilesOnSWD = [];
                                argsAsObject = lodash_1.default.fromPairs(stepToRun.args.map(function (x) { return [x.injectInto, x.value]; }));
                                fs_extra_1.outputJsonSync(path_1.resolve(swdOnHost_1, config_1.default.argsFileOnSWD), argsAsObject);
                                taskRelatedFilesOnSWD.push({ path: config_1.default.argsFileOnSWD, opt: 'ro' });
                                return [4 /*yield*/, knex('uploaded_file')
                                        .select(['id', 'inject_into', 'meta'])
                                        .where('step_id', stepToRun.id)];
                            case 9:
                                currentUploadedFiles = _g.sent();
                                stepUploadedFilesDir = path_1.resolve(config_1.default.stepPathBaseOnHost, stepToRun.id, 'uploaded-files');
                                fs_extra_1.ensureDirSync(stepUploadedFilesDir);
                                _b = (_a = taskRelatedFilesOnSWD.push).apply;
                                _c = [taskRelatedFilesOnSWD];
                                return [4 /*yield*/, Promise.all(currentUploadedFiles.map(function (upl) { return __awaiter(_this, void 0, void 0, function () {
                                        var fpGbox;
                                        return __generator(this, function (_a) {
                                            fpGbox = path_1.resolve(swdOnHost_1, config_1.default.uploadedFilesDirOnSWD, upl.inject_into, upl.meta.name);
                                            console.log(chalk_1.default.bold.blue("linking uploadedFile:  " + path_1.resolve(config_1.default.uploadPathOnHost, upl.id) + "  ->  " + fpGbox));
                                            fs_extra_1.removeSync(fpGbox);
                                            fs_extra_1.ensureLinkSync(path_1.resolve(config_1.default.uploadPathOnHost, upl.id), fpGbox);
                                            return [2 /*return*/, { path: path_1.join(config_1.default.uploadedFilesDirOnSWD, upl.inject_into, upl.meta.name), opt: 'ro' }];
                                        });
                                    }); }))];
                            case 10:
                                _b.apply(_a, _c.concat([(_g.sent())]));
                                return [4 /*yield*/, knex('import')
                                        .select(['id', 'export_id', 'inject_into'])
                                        .where('step_id', stepToRun.id)];
                            case 11:
                                currentImports = _g.sent();
                                _e = (_d = taskRelatedFilesOnSWD.push).apply;
                                _f = [taskRelatedFilesOnSWD];
                                return [4 /*yield*/, Promise.all(currentImports.map(function (imp) { return __awaiter(_this, void 0, void 0, function () {
                                        var fpGbox;
                                        return __generator(this, function (_a) {
                                            fpGbox = path_1.resolve(swdOnHost_1, config_1.default.importsDirOnSWD, imp.inject_into);
                                            console.log(chalk_1.default.bold.blue("linking import:  " + path_1.resolve(config_1.default.dataPathOnHost, imp.export_id) + "  ->  " + fpGbox));
                                            fs_extra_1.removeSync(fpGbox);
                                            fs_extra_1.ensureLinkSync(path_1.resolve(config_1.default.dataPathOnHost, imp.export_id), fpGbox);
                                            return [2 /*return*/, { path: path_1.join(config_1.default.importsDirOnSWD, imp.inject_into), opt: 'ro' }];
                                        });
                                    }); }))];
                            case 12:
                                _e.apply(_d, _f.concat([(_g.sent())]));
                                exportsDirOnHost_1 = path_1.resolve(swdOnHost_1, config_1.default.exportsDirOnSWD);
                                fs_extra_1.ensureDirSync(exportsDirOnHost_1);
                                taskRelatedFilesOnSWD.push({ path: config_1.default.exportsDirOnSWD, opt: '' });
                                exportsAnnoFileOnHost = path_1.resolve(swdOnHost_1, config_1.default.exportsAnnoFileOnSWD);
                                fs_extra_1.ensureFileSync(exportsAnnoFileOnHost);
                                taskRelatedFilesOnSWD.push({ path: config_1.default.exportsAnnoFileOnSWD, opt: '' });
                                resultsFileOnHost_1 = path_1.resolve(swdOnHost_1, config_1.default.resultsFileOnSWD);
                                fs_extra_1.ensureFileSync(resultsFileOnHost_1);
                                taskRelatedFilesOnSWD.push({ path: config_1.default.resultsFileOnSWD, opt: '' });
                                debugDirOnHost = path_1.resolve(swdOnHost_1, config_1.default.debugDirOnSWD);
                                fs_extra_1.ensureDirSync(debugDirOnHost);
                                taskRelatedFilesOnSWD.push({ path: config_1.default.debugDirOnSWD, opt: '' });
                                return [4 /*yield*/, getExeCmd_1.getExeCmdArgv(knex, stepToRun.gbox, taskRelatedFilesOnSWD, swdOnHost_1)];
                            case 13:
                                exeCmdArgv = _g.sent();
                                gboxStdoutFileOnHost = path_1.resolve(swdOnHost_1, 'stdout');
                                fs_extra_1.ensureFileSync(gboxStdoutFileOnHost);
                                gboxStderrFileOnHost = path_1.resolve(swdOnHost_1, 'stderr');
                                fs_extra_1.ensureFileSync(gboxStderrFileOnHost);
                                gboxStartTime = process.hrtime();
                                return [4 /*yield*/, runGbox_1.default(knex, exeCmdArgv, gboxStdoutFileOnHost, gboxStderrFileOnHost, function () { return __awaiter(_this, void 0, void 0, function () {
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0: return [4 /*yield*/, knex('step')
                                                        .select(['status'])
                                                        .where({ id: stepToRun.id })];
                                                case 1: return [2 /*return*/, (_a.sent())[0].status === 'interception_requested'];
                                            }
                                        });
                                    }); }, function (obj) { return __awaiter(_this, void 0, void 0, function () {
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0: return [4 /*yield*/, knex('step')
                                                        .update({ 'errors': obj })
                                                        .where({ id: stepToRun.id })];
                                                case 1:
                                                    _a.sent();
                                                    return [2 /*return*/];
                                            }
                                        });
                                    }); })];
                            case 14:
                                _g.sent();
                                return [4 /*yield*/, knex('step')
                                        .update({ 'errors': null })
                                        .where({ id: stepToRun.id })];
                            case 15:
                                _g.sent();
                                gboxRunningDuration_1 = process.hrtime(gboxStartTime);
                                console.log("Gbox has finished running in " + gboxRunningDuration_1[0] + " seconds and " + gboxRunningDuration_1[1] + " nanoseconds");
                                results_1 = (function () {
                                    try {
                                        return fs_extra_1.readJsonSync(resultsFileOnHost_1);
                                    }
                                    catch (e) {
                                        throw new Error("Error reading " + resultsFileOnHost_1 + ".");
                                    }
                                })();
                                annotatedExportsFromJsonFile = (function () {
                                    try {
                                        return fs_extra_1.readJsonSync(path_1.resolve(swdOnHost_1, config_1.default.exportsAnnoFileOnSWD));
                                    }
                                    catch (e) {
                                        console.log("No dynamic export annotation: " + config_1.default.exportsAnnoFileOnSWD + " does not exist or is not parsed correctly.");
                                        return [];
                                    }
                                })();
                                return [4 /*yield*/, Promise.all(annotatedExportsFromJsonFile.map(function (exp) { return __awaiter(_this, void 0, void 0, function () {
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0:
                                                    if (exp.extractFrom == null) {
                                                        throw new Error("The field \"extractFrom\" not found. This field is required in a dynamic export annotation.");
                                                    }
                                                    return [4 /*yield*/, knex('export').insert({
                                                            step_id: stepToRun.id,
                                                            extract_from: exp.extractFrom,
                                                            kind: exp.kind || null,
                                                            meta: exp.meta || null,
                                                        })];
                                                case 1:
                                                    _a.sent();
                                                    return [2 /*return*/];
                                            }
                                        });
                                    }); }))];
                            case 16:
                                _g.sent();
                                return [4 /*yield*/, knex('export')
                                        .select(['id', 'extract_from'])
                                        .where({ step_id: stepToRun.id })];
                            case 17:
                                annotatedExportsFromDatabase = _g.sent();
                                return [4 /*yield*/, Promise.all(annotatedExportsFromDatabase.map(function (exp) { return __awaiter(_this, void 0, void 0, function () {
                                        var expFile;
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0:
                                                    expFile = path_1.resolve(exportsDirOnHost_1, exp.extract_from);
                                                    if (!fs_extra_1.pathExistsSync(expFile)) {
                                                        throw new Error("The annotated export file " + expFile + " does not exist.");
                                                    }
                                                    fs_extra_1.removeSync(path_1.resolve(config_1.default.dataPathOnHost, exp.id));
                                                    console.log(chalk_1.default.bold.blue("linking export:  " + expFile + "  ->  " + path_1.resolve(config_1.default.dataPathOnHost, exp.id)));
                                                    fs_extra_1.ensureLinkSync(expFile, path_1.resolve(config_1.default.dataPathOnHost, exp.id));
                                                    return [4 /*yield*/, knex('export')
                                                            .update({ is_populated: true })
                                                            .where({ id: exp.id })];
                                                case 1:
                                                    _a.sent();
                                                    return [2 /*return*/];
                                            }
                                        });
                                    }); }))];
                            case 18:
                                _g.sent();
                                /**
                                 * Finally we writes results and errors into the database
                                 */
                                return [4 /*yield*/, knex.transaction(function (trx) { return __awaiter(_this, void 0, void 0, function () {
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0: return [4 /*yield*/, trx('step')
                                                        .update({
                                                        status: 'done',
                                                        results: JSON.stringify(results_1),
                                                        state: JSON.stringify({
                                                            executionDuration: { second: gboxRunningDuration_1[0], nanoSecond: gboxRunningDuration_1[1] },
                                                        }),
                                                    })
                                                        .where('id', stepToRun.id)];
                                                case 1:
                                                    _a.sent();
                                                    return [2 /*return*/];
                                            }
                                        });
                                    }); })];
                            case 19:
                                /**
                                 * Finally we writes results and errors into the database
                                 */
                                _g.sent();
                                return [3 /*break*/, 22];
                            case 20:
                                e_2 = _g.sent();
                                swdOnHost = path_1.resolve(config_1.default.stepPathBaseOnHost, stepToRun.id);
                                gboxStdout = utils_1.safeReadFileSync(path_1.resolve(swdOnHost, 'stdout')).toString();
                                gboxStderr = utils_1.safeReadFileSync(path_1.resolve(swdOnHost, 'stderr')).toString();
                                gboxLogMessage_1 = [
                                    gboxStdout ? "---stdout---\n" + gboxStdout : '---stdout is empty---',
                                    gboxStderr ? "---stderr---\n" + gboxStderr : '---stderr is empty---',
                                ].join('\n');
                                return [4 /*yield*/, knex.transaction(function (trx) { return __awaiter(_this, void 0, void 0, function () {
                                        return __generator(this, function (_a) {
                                            switch (_a.label) {
                                                case 0: 
                                                // TODO: change to reset_step force
                                                return [4 /*yield*/, trx.raw('SELECT reset_step_recursively(?)', [stepToRun.id])];
                                                case 1:
                                                    // TODO: change to reset_step force
                                                    _a.sent();
                                                    // noinspection JSReferencingMutableVariableFromClosure
                                                    return [4 /*yield*/, trx('step')
                                                            .update({
                                                            results: null,
                                                            errors: JSON.stringify(config_1.default.__DEV__
                                                                ? [{ source: 'NodeJS', message: e_2.stack }, { source: 'Gbox', message: gboxLogMessage_1 }]
                                                                : [{ source: 'Granatum', messsage: 'An internal error just happened. Please try again.' }]),
                                                        })
                                                            .where('id', stepToRun.id)];
                                                case 2:
                                                    // noinspection JSReferencingMutableVariableFromClosure
                                                    _a.sent();
                                                    return [2 /*return*/];
                                            }
                                        });
                                    }); })];
                            case 21:
                                _g.sent();
                                console.error(e_2);
                                return [3 /*break*/, 22];
                            case 22: return [2 /*return*/];
                        }
                    });
                };
                _a.label = 1;
            case 1:
                if (!true) return [3 /*break*/, 3];
                return [5 /*yield**/, _loop_1()];
            case 2:
                _a.sent();
                return [3 /*break*/, 1];
            case 3: return [2 /*return*/];
        }
    });
}); });
