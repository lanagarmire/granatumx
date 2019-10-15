export default {
  port: +process.env.PORT || 34567,
  databaseUrl: process.env.DATABASE_URL || 'postgres://postgres:12qw@localhost:5433/granatum',
  uploadPath: '/var/granatum/uploaded_files',
  dataPath: '/var/granatum/data',
  downloadFilePath: '/var/granatum/download_files',

  // api: {
  //   clientUrl: process.env.API_CLIENT_URL || '',
  //   serverUrl: process.env.API_SERVER_URL || `http://localhost:${process.env.PORT || 3000}`,
  // },
  // analytics: {
  //   googleTrackingId: process.env.GOOGLE_TRACKING_ID, // UA-XXXXX-X
  // },
  auth: {
    jwt: { secret: process.env.JWT_SECRET || 'ufdMnOWZ41' },
    // facebook: {
    //   id: process.env.FACEBOOK_APP_ID || '186244551745631',
    //   secret: process.env.FACEBOOK_APP_SECRET || 'a970ae3240ab4b9b8aae0f9f0661c6fc',
    // },
    // google: {
    //   id: process.env.GOOGLE_CLIENT_ID || '251410730550-ahcg0ou5mgfhl8hlui1urru7jn5s12km.apps.googleusercontent.com',
    //   secret: process.env.GOOGLE_CLIENT_SECRET || 'Y8yR9yZAhm9jQ8FKAL8QIEcd',
    // },
    // twitter: {
    //   key: process.env.TWITTER_CONSUMER_KEY || 'Ie20AZvLJI2lQD5Dsgxgjauns',
    //   secret: process.env.TWITTER_CONSUMER_SECRET || 'KTZ6cxoKnEakQCeSpZlaUCJWGAlTEBJj0y2EMkUBujA7zWSvaQ',
    // },
  },
};
