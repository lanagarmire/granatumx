import fs from 'fs';
import https from 'https';

export default function startHttpsServer(httpApp, cert, key) {
  const options = {
    cert: fs.readFileSync(cert),
    key: fs.readFileSync(key),
  };

  return https.createServer(options, httpApp).listen(443);
}
