# Application hubs guide

A note on the terminology: We refer to singleuser as being the software stack that the user directly interacts with. The _singleusers_ in _astrohub_ are different from the singular _singleuser_ used in `cyberhubs/corehub`. The _astrohub_ _singleuser_ images use the `cyberhubs/corehub` singleuser image to build. This means that before the _astrohub_ _singleusers_ are deployed, the corehub _singleuser_ and _multiuser_ must be built. See GitHub `cyberlaboratories/cyberhubs` for access to corehub.  

## User Guide
To deploy the singleuser image of the desired _astrohub_, navigate to the directory of the singleuser you wish to employ (_wendihub_, _mlhub_, _mesahub_, etc).

Next step is to configure the `jupyterhub-config-script.sh`, specifically, the environment variable `JUPYTER_SGLEUSR_IMG` should be set to the name of the singleuser you would like. Typically, choose a name like `cyberhubs/wendihub` or `cyberhubs/nameofsingleuser`. To set environment variables once script is edited, just use the 
```source jupyterhub-config-script.sh``` 
command.

Now deploy the image by entering 
```make build``` 
and the image should build with the given name in `JUPYTER_SGLEUSR_IMG`. 

## Developer Guide

### Trusted Notebooks

By default, the notebooks developed in the _astrohubs_ _singleuser_ images are untrusted by jupyterhub as a security measure. This means that every time a new notebook is added, its path must be added to `trusted_nb.sh`.

Thus, if you are developing notebooks that you would like to permanently exist in the _singleuser_ environment, you must make this change in order for the notebooks to automatically run. If you are adding notebooks permanently to the _singleuser_ image you have built, this implies you have made changes to the `Dockerfile` to load the notebook into the _singleuser_ environment, and a rebuild of the `Dockerfile` is necessary following the next step.   

Each line in `trusted_nb.sh` is a separate notebook in the _singluser_ environment. The command is simply 
```jupyter trust /path/to/notebook```
Then a rebuild of the _singleuser_ image is necessary in order for these changes to be made.  