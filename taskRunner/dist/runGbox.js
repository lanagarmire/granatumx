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
var child_process_1 = require("child_process");
var fs_extra_1 = require("fs-extra");
var config_1 = __importDefault(require("./config"));
exports.default = (function (knex, exeCmdArgv, gboxStdoutFileOnHost, gboxStderrFileOnHost, checkKillReq, writeMessageToErrors) {
    return new Promise(function (rsv, rjt) {
        var gboxProcess = child_process_1.spawn(exeCmdArgv[0], exeCmdArgv.slice(1), { stdio: 'pipe' });
        gboxProcess.stdout.pipe(fs_extra_1.createWriteStream(gboxStdoutFileOnHost));
        gboxProcess.stderr.pipe(fs_extra_1.createWriteStream(gboxStderrFileOnHost));
        var stdoutStr = '';
        gboxProcess.stdout.on('data', function (chunk) { return __awaiter(_this, void 0, void 0, function () {
            return __generator(this, function (_a) {
                switch (_a.label) {
                    case 0:
                        stdoutStr += chunk;
                        stdoutStr = stdoutStr.slice(-10000);
                        if (!(writeMessageToErrors != null)) return [3 /*break*/, 2];
                        return [4 /*yield*/, writeMessageToErrors({ stdout: stdoutStr, stderr: stderrStr })];
                    case 1:
                        _a.sent();
                        _a.label = 2;
                    case 2: return [2 /*return*/];
                }
            });
        }); });
        var stderrStr = '';
        gboxProcess.stderr.on('data', function (chunk) { return __awaiter(_this, void 0, void 0, function () {
            return __generator(this, function (_a) {
                switch (_a.label) {
                    case 0:
                        stderrStr += chunk;
                        stderrStr = stderrStr.slice(-10000);
                        if (!(writeMessageToErrors != null)) return [3 /*break*/, 2];
                        return [4 /*yield*/, writeMessageToErrors({ stdout: stdoutStr, stderr: stderrStr })];
                    case 1:
                        _a.sent();
                        _a.label = 2;
                    case 2: return [2 /*return*/];
                }
            });
        }); });
        console.log("Gbox has started. (Timelimit = " + config_1.default.GBOX_TIME_LIMIT + " ms)");
        // kill the process if it's taken too long
        var killTimeoutId = setTimeout(function () { return __awaiter(_this, void 0, void 0, function () {
            return __generator(this, function (_a) {
                gboxProcess.kill('SIGKILL');
                return [2 /*return*/];
            });
        }); }, +config_1.default.GBOX_TIME_LIMIT);
        // kill the process if requested by the user
        var killInterceptedId = setInterval(function () { return __awaiter(_this, void 0, void 0, function () {
            var _a;
            return __generator(this, function (_b) {
                switch (_b.label) {
                    case 0:
                        _a = checkKillReq;
                        if (!_a) return [3 /*break*/, 2];
                        return [4 /*yield*/, checkKillReq()];
                    case 1:
                        _a = (_b.sent());
                        _b.label = 2;
                    case 2:
                        if (_a) {
                            gboxProcess.kill('SIGKILL');
                        }
                        return [2 /*return*/];
                }
            });
        }); }, 1000);
        gboxProcess.on('exit', function (code, signal) {
            clearInterval(killInterceptedId);
            clearTimeout(killTimeoutId);
            if (code === 0) {
                rsv();
            }
            else {
                rjt(new Error("The gbox process exited with the code " + code + ", caused by signal " + JSON.stringify(signal)));
            }
        });
    });
});
