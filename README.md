
![logo](https://raw.githubusercontent.com/Novartis/railroadtracks/master/doc/_static/logo_rrt.png)

Railroadtracks
==============

Railroadtracks is a Python package to handle connected computation steps for DNA and RNA Seq.

[![PyPI version](https://badge.fury.io/py/railroadtracks.svg)](http://badge.fury.io/py/railroadtracks)

| Branch        |                                                                                                                                                   |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| Master        | [![Build Status](https://travis-ci.org/Novartis/railroadtracks.svg?branch=master)](https://travis-ci.org/Novartis/railroadtracks/branches)        |
| version_0.3.x | [![Build Status](https://travis-ci.org/Novartis/railroadtracks.svg?branch=version_0.3.x)](https://travis-ci.org/Novartis/railroadtracks/branches) |
| version_0.4.x | [![Build Status](https://travis-ci.org/Novartis/railroadtracks.svg?branch=version_0.4.x)](https://travis-ci.org/Novartis/railroadtracks/branches) |



Installation
------------

The documentation should be consulted for information about the requirements and the installation process.

While we are working on a link to a build of the documentation, the Sphinx source for it
are be consulted: [doc/installation.rst](https://raw.githubusercontent.com/Novartis/railroadtracks/master/doc/installation.rst).

Released versiond are available on Pypi. Installing the latest release is as easy as:
```
pip install railroadtracks
```

The use the development version, the master branch can be installed with
```
pip install https://github.com/Novartis/railroadtracks/archive/master.zip
```

A [Docker](http://www.docker.io) container is also availale (`lgautier/railroadtracks`).

```
# pull the Docker images from docker.io
docker pull lgautier/railroadtracks 
```

```
# Run an ipython notebook server from the container
docker run -it -p 8888:8888 -w /usr/local/packages/railroadtracks/doc/notebooks lgautier/railroadtracks ipython notebook --ip=0.0.0.0 --no-browser
```

A browser running on the same machine as the container can then be pointed to:
[http://localhost:8888](http://localhost:8888)

If on a non-Linux system using [boot2docker](http://boot2docker.io), you may have to point your browser to the IP mapped to the VM docker is running in. The IP can be obtained with:

```
boot2docker ip
```

Documentation
-------------

Automated building of the documentation is being worked on.

In the meantime there is a [tutorial as an ipython notebook](http://nbviewer.ipython.org/github/Novartis/railroadtracks/blob/master/doc/notebooks/railroadtracks_tutorial.ipynb")


