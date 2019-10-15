# GranatumX

<!-- The next generation of the most popular graphical single-cell RNA-seq analysis pipeline for genomics scientists. -->

Welcome to the next generation of graphical single-cell RNA-seq analysis pipeline for genomics scientists!

## What's new

<!--
The entire architecture has been re-designed from the ground up. Here are some of the killing
features of GranatumX:

  * A persistent user management system that securely stores all your data
  * A project-system where you can have multiple projects with different configurations
  * A complete de-coupled front-end and back-end. Deploy it everywhere!
  * A modern plugin-system. Put your favorite library in Granatum today!
  * New interface! Clean and intuitive as always, but with much more functionalities!
-->

* Customize processing with a modern plug-in system – add your favorite library!
* Leverage multiple servers by deploying a single backend database with multiple frontend clients for processing
* Save and manage multiple projects from one secure personal account
* New clean and intuitive interface with much more functionality

## Architecture

<!-- Web interface: -->

The entire architecture has been re-designed from the ground up

* Server-side:
  * NodeJS + ExpressJS (web-serving)
  * PostgreSQL + Knex (database)
  * Graphile (postgres to graphql schema)
  * React (server-side rendering)
  * Docker (containerized execution)
* Client-side:
  * Apollo (graphql client)
  * React (component rendering)
  * Material-UI (theme)

<!-- 
      - PostgreSQL (database)
-->

<!--
Task scheduling and executing:
  * NodeJS + Knex + PostgreSQL
  * Docker
-->

## Creating a Gbox

GranatumX ships with a template to make it easier for developers to create new Gboxes. Below is a checklist for creating a new Gbox.

1. Download the GranatumX Source Code `git clone https://gitlab.com/xz/GranatumX`
2. In the GranatumX folder, run `make setup`. After setup completes, GranatumX should start on port `34567`
3. Copy the Gbox template into the `g_packages` directory: `cp gboxTemplate g_packages/yourPackageName`
4. Edit your `package.yaml` file with your package information. See other packages in the `g_packages` folder for more examples.
5. Edit your `Dockerfile` by adding any package installation scripts your Gbox requires.
6. Edit your `main.py` file with your application code (or whatever file you specified in `package.yaml > gboxes > endpoints > backend > cmd`). The Python package `granatum_sdk` contains helper methods for easily interacting with the GranatumX core.
7. Install your new gbox with `cd gboxes; yarn run installEverything`
8. Refresh GranatumX and test your new Gbox.

## Running your own server

### Environment Variables

If starting the server with sudo you may need to add the -E flag to use your current environment variables (sudo -E make start).

* PORT: Port # to listen to, i.e. 80, 3000, 8888
* DATABASE_URL: Url to connect to your Postgresql database server
* SSL_CERT: (optional) To serve Granatum on port 443 over SSL, set to certificate filepath
* SSL_KEY: (optional) To serve Granatum on port 443 over SSL, set to private key filepath


