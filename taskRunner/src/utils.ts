import { readFileSync } from 'fs';
import { resolve } from 'path';

export const timeout = (t: number) =>
  new Promise((rsv) => {
    setTimeout(() => {
      rsv();
    }, t);
  });

export const safeReadFileSync = (path?: string) => {
  if (path == null) {
    return new Buffer(0);
  }

  try {
    return readFileSync(path);
  } catch (e) {
    return new Buffer(0);
  }
};
