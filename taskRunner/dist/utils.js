"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var fs_1 = require("fs");
exports.timeout = function (t) {
    return new Promise(function (rsv) {
        setTimeout(function () {
            rsv();
        }, t);
    });
};
exports.safeReadFileSync = function (path) {
    if (path == null) {
        return new Buffer(0);
    }
    try {
        return fs_1.readFileSync(path);
    }
    catch (e) {
        return new Buffer(0);
    }
};
