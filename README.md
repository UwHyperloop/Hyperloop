<img src="http://www.washington.edu/brand/files/2014/09/Signature_Left_Purple_Hex.png" height="30" />

Updated version of the Hyperloop OpenMDAO plugin, written for OpenMDAO 1.0.

# Hyperloop Model


This is an open-source system model of the Hyperloop transportation 
system outlined in Elon Musk’s proposal. It's primary 
purpose is to provide a completely open-source multidisciplinary model 
as a central point for continued crowd-sourced refinement of the concept. 
The current model focuses mainly on the compression system: inlet compressors, 
heat exchangers, ducting, and nozzle. The entire model is a work in progress 
and was built with the hope that external participation can help both 
refine our estimates and expand the model to include other important 
aspects of the hyperloop concept. We've put together [detailed documentation](http://openmdao-plugins.github.io/Hyperloop/index.html) of 
the models and some conclusions we've drawn with them. 

![Hyperloop Rendering](https://raw.github.com/OpenMDAO-Plugins/Hyperloop/master/docs/images/hyperloop.png "Hyperloop Rendering with a design speed of Mach .8")

Due to the cross-disciplinary nature of the problem, an overarching framework is 
needed to orchestrate the interaction between models of the various subsystems. So 
the entire model is built ontop of the OpenMDAO framework. 
OpenMDAO is an open source effort, out of the NASA Glenn Research Center,  
that suits the needs of this problem well. It provides advanced modeling and optimization 
capabilities and provides a flexible structure to grow and refine the hyperloop models. 
OpenMDAO includes a plugin system, which we used to build this model. The entire thing is built 
and distributed as an OpenMDAO plugin, To further improve the ability for others to expand and improve on this work. 

Our initial results indicate that, while the concept is still very viable, the tube size will need 
to be around twice as large as originally proposed by the Tesla/SpaceX team. The tube diameter 
is closely tied to the maximum speed the pod travels at. 

![Tube Diameter vs Mach](https://github.com/OpenMDAO-Plugins/Hyperloop/blob/master/docs/images/mach_vs_rad.png?raw=true "Hyperloop Tube Diameter vs Mach Number")


## Pre-Reqs

This is an [OpenMDAO Plugin](http://openmdao.org/) so you need to install OpenMDAO 1.2.0 or higher, PyCycle2 and all of their dependencies to use it.

### OpenMDAO

If you need to install OpenMDAO follow the [install instructions](http://openmdao.org/docs/getting-started/index.html). 

If you're using a Mac, then your best bet for installing all the pre-requisites is to use 
homebrew. You can follow these [detailed instructions](http://www.lowindata.com/2013/installing-scientific-python-on-mac-os-x/)
but once you have homebrew installed and setup, here is the short version: 

```
brew install git
brew install python
brew install gfortran
pip install numpy
pip install scipy
brew install freetype
pip install matplotlib
```

### PyCycle
The hyperloop model depends on another plugin, [PyCycle](https://github.com/OpenMDAO-Plugins/pyCycle) a thermodynamic cycle modeling tool.
You should follow the installation instruction for that plugin before moving on.  


## Installation
Then, make sure you have activated your OpenMDAO environment. Now you should clone this 
repository (or your fork of it) to your local machine.  Either way. 

    git clone https://github.com/OpenMDAO-Plugins/Hyperloop.git

Next navigate to the top level of the repo and issue the following command 

    python setup.py develop

This will install plugin and let you make modifications as you like. 


## Please Read The Docs
<s>You can read the [online version of the docs](http://openmdao-plugins.github.io/Hyperloop/), which tracks the latest version of the code on https://github.com/OpenMDAO-Plugins/Hyperloop. Or you can read the docs directly from the repository. Once you download the repo, you need build the docs one time. If you make any changes to them, you need to rebuild them to see the updates. To build the docs navigate to the top of the repository and use the following command to build the docs.
    
    plugin build docs

Then, from anywhere in an activated environment you can open the docs with the following command

    plugin docs hyperloop

The docs give a lot of background on the model, and explain the thinking that went into its 
structure. They also give some usage examples.</s> If you plan to make any of your own contributions 
to this work, you might find the information in the [OpenMDAO developer docs](http://openmdao.org/docs/dev-guide/index.html) 
useful.

### Quickstart
If you just can't wait to get started, the simulation can be run straight out of the box with:

    cd src/hyperloop
    python hyperloop_sim.py


### Credit
These models were created by Justin Gray and Jeffrey Chin, with contributions from Scott Jones, Jeff Berton, and Chris Heath. The models have been expanded and modified by Brent Schroeter and the UWashington Hyperloop team.

### License
This model is covered by an Apache V2.0 License




