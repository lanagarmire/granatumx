import { spawn } from 'child_process';
import { createWriteStream } from 'fs-extra';
import * as Knex from 'knex';
import config from './config';

export default (
  knex: Knex,
  exeCmdArgv: string[],
  gboxStdoutFileOnHost: string,
  gboxStderrFileOnHost: string,
  checkKillReq?: () => Promise<boolean>,
  writeMessageToErrors?: (obj: object) => Promise<void>,
) =>
  new Promise((rsv, rjt) => {
    const gboxProcess = spawn(exeCmdArgv[0], exeCmdArgv.slice(1), { stdio: 'pipe' });
    gboxProcess.stdout.pipe(createWriteStream(gboxStdoutFileOnHost));
    gboxProcess.stderr.pipe(createWriteStream(gboxStderrFileOnHost));

    let stdoutStr = '';
    gboxProcess.stdout.on('data', async (chunk) => {
      stdoutStr += chunk;
      stdoutStr = stdoutStr.slice(-10000);
      if (writeMessageToErrors != null) {
        await writeMessageToErrors({ stdout: stdoutStr, stderr: stderrStr });
      }
    });

    let stderrStr = '';
    gboxProcess.stderr.on('data', async (chunk) => {
      stderrStr += chunk;
      stderrStr = stderrStr.slice(-10000);
      if (writeMessageToErrors != null) {
        await writeMessageToErrors({ stdout: stdoutStr, stderr: stderrStr });
      }
    });

    console.log(`Gbox has started. (Timelimit = ${config.GBOX_TIME_LIMIT} ms)`);

    // kill the process if it's taken too long
    const killTimeoutId = setTimeout(async () => {
      gboxProcess.kill('SIGKILL');
    }, +config.GBOX_TIME_LIMIT);

    // kill the process if requested by the user
    const killInterceptedId = setInterval(async () => {
      if (checkKillReq && (await checkKillReq())) {
        gboxProcess.kill('SIGKILL');
      }
    }, 1000);

    gboxProcess.on('exit', (code, signal) => {
      clearInterval(killInterceptedId);
      clearTimeout(killTimeoutId);
      if (code === 0) {
        rsv();
      } else {
        rjt(new Error(`The gbox process exited with the code ${code}, caused by signal ${JSON.stringify(signal)}`));
      }
    });
  });
