/* eslint-disable no-console */
// Run with: npx babel-node status.js

import chalk from 'chalk';
import { spawn, spawnSync } from 'child_process';
import Docker from 'dockerode';
import find from 'find-process';
import fs from 'fs';
import glob from 'glob';
import yaml from 'js-yaml';
import fetch from 'node-fetch';

import Knex from './src/server/getSingletonKnex';

const REFRESH_INTERVAL = 5000;

// Todo: add status of all the gboxes

const isDockerUp = async () => {
  try {
    const docker = new Docker({ socketPath: '/var/run/docker.sock' });
    await docker.listContainers();
    return true;
  } catch (err) {
    return false;
  }
};

const isPostgresUp = async (knex) => {
  try {
    // Run arbitrary query to check status.
    await knex.select('*').from('gbox');
    return true;
  } catch (err) {
    return false;
  }
};

const isProcessRunning = async (name) => {
  try {
    const list = await find('name', name);
    return list.length > 0;
  } catch (err) {
    return false;
  }
};

const getGboxName = (filename) => filename.replace(/^.*\/([^/]+)\.yaml$/, '$1');

const getDockerCommandToRun = (filename) => {
  const spec = yaml.safeLoad(fs.readFileSync(filename, 'utf8')) as any;
  const backend = spec.endpoints.backend;
  const command = [
    `docker run -i --rm --name ${getGboxName(filename)} -e GRANATUM_SWD=/data`,
    backend.image,
    backend.cmd,
  ].join(' ');
  console.log(command);
  return command;
};

// todo: debug. why aren't some docker instances closing?
// https://hackernoon.com/another-reason-why-your-docker-containers-may-be-slow-d37207dec27f
const testGboxAsync = (filename) =>
  new Promise((resolve) => {
    try {
      const command = getDockerCommandToRun(filename);
      const commandWords = command.split(' ');
      const execution = spawn(commandWords[0], commandWords.slice(1));

      execution.stderr.on('data', (err) => {
        console.log(`${getGboxName(filename)} stderr: ${err}`);
      });

      execution.on('close', (exitCode) => {
        if (exitCode > 0) {
          console.log(`${getGboxName(filename)} exitCode ${exitCode}`);
          resolve(false);
        }
        resolve(true);
      });
    } catch (err) {
      console.log(err);
      resolve(false);
    }
  });

const testGbox = async (filename) => {
  try {
    const command = getDockerCommandToRun(filename);
    const commandWords = command.split(' ');
    const execution = spawnSync(commandWords[0], commandWords.slice(1));
    const exitCode = execution.status;

    if (exitCode > 0) {
      throw new Error(execution.stderr.toString());
    }

    return true;
  } catch (err) {
    //  console.log(err);
    return false;
  }
};

const printTest = (testName, testResult, suffix = `up`) => {
  console.log(testResult ? chalk.green(`${testName} is ${suffix}`) : chalk.red(`${testName} is NOT ${suffix}`));
};

const testGboxes = async () => {
  const gboxYamls = glob.sync('../gboxes/gboxSpecs/*.yaml');
  const promises = gboxYamls.map(async (filename) =>
    testGboxAsync(filename).then((result) => {
      printTest(`Gbox: ${getGboxName(filename)}`, result);
    }),
  );

  await Promise.all(promises);
};

const runTests = async () => {
  const knex = await Knex();
  printTest('Postgres/Node Connection', await isPostgresUp(knex));
  printTest('Docker', await isDockerUp());
  printTest('TaskRunner', await isProcessRunning('taskRunner'));
  printTest('WebApp', await isProcessRunning('granatumWebApp'));

  // Todo: be smarter about what ports to check
  const ports = [80, 443, 3000, 34567];

  ports.forEach(async (port) => {
    try {
      const url = port === 443 ? 'https://granatum.garmiregroup.org' : `http://localhost:${port}/index.html`;
      const res = await fetch(url, {
        redirect: 'manual',
      });
      printTest('WebApp', true, `listening on port ${port}`);
      const text = await res.text();
      printTest('WebApp', text.includes('Granatum'), `working properly on port ${port}`);
    } catch (err) {
      printTest('WebApp', false, `listening on port ${port}`);
    }
  });

  if (process.argv.includes('--gbox')) {
    testGboxes();
  }

  knex.destroy();
};

const delay = (t) =>
  new Promise((resolve) => {
    setTimeout(() => {
      resolve();
    }, t);
  });

setImmediate(async () => {
  const animateFrames = '\\-/|'.split('');
  let animateFrameNum = 0;
  // eslint-disable-next-line no-constant-condition
  while (true) {
    // eslint-disable-next-line
    console.log('\x1B[2J\x1B[0f');
    animateFrameNum = (animateFrameNum + 1) % animateFrames.length;
    console.log(`${animateFrames[animateFrameNum]} Refreshing in ${REFRESH_INTERVAL / 1000} seconds ...`);
    console.log('');

    runTests();

    // eslint-disable-next-line no-await-in-loop
    await delay(REFRESH_INTERVAL);
  }
});
