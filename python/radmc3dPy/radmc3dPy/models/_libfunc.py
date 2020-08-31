from __future__ import absolute_import
from __future__ import print_function

import subprocess as sp
import os
import importlib
import warnings

def getTemplateModel():
    """Create a copy of the template model file in the current working directory. 

    The PYTHONPATH environment variable is checked for the installation path of radmc3dPy and
    the template file is copied from the first hit in the path list. 
    """
   
    # Get the installation directory of the radmc3dPy package
    import radmc3dPy
    mod_path = radmc3dPy.__file__.strip()[:-12]

    # Copy the model template to the current working directory
    command = 'cp -v '+mod_path+'/template.py '+os.getcwd()
    try:
        os.system(command)
    except Exception as ex:
        print(ex)

        return False
    return True


def updateModelList():
    """Updates the list of models in the library.

    This function gets the names of all .py files in the radmc3dPy/models directory. Then tries to import them as
    a rought check for consistency. If the content of the file can be imported, the name of the model file is added
    to the model list. The model modules are then added to the _modellist.py in the radmc3dPy/models directory and
    can be used from that point on. 

    NOTE 1 : Write permissions is needed to update the list of models as the radmc3dPy/models/_modellist.py file
        has to be re-written.

    NOTE 2 : After the updateModelList has been run the module must be either reloaded or one has to exit the current
        python session and relaunch it to see the effect of the updated model list. 
    """
        
    mod_names = []
   
    # Get the name of all model files in the module directory
    import radmc3dPy
    mod_path = radmc3dPy.__file__.strip()[:-12]
    dum = sp.Popen(['ls -1 '+mod_path+'/models/*.py'], shell=True,
                   stdout=sp.PIPE, stderr=sp.PIPE).communicate()[0].split()
    homeDir = os.getcwd()

    os.chdir(mod_path+'/models')
    for i in range(len(dum)):
        sdum = dum[i].split("/")
        modname = sdum[len(sdum)-1][:-3]

        if (modname!='template')&(modname!='__init__')&(modname!='_libfunc')&(modname!='_modellist'):
            # Try to import the model. If it can be imported without any error/exception add the model
            #  to the list of models
           
            mod_names.append(modname)
            try:
                mdl = importlib.import_module(modname)
            except ImportError:
                try:
                    mdl = importlib.import_module('radmc3dPy.models.'+modname)
                except ImportError:
                    warnings.warn(modname + ' cannot be imported. Skipping ...', RuntimeWarning)
                    mod_names.remove(modname)
    

    os.chdir(homeDir)
    
    with open(mod_path+'/models/_modellist.py', 'w') as wfile:

        for imod in range(len(mod_names)):
                wfile.write('import %s\n'%mod_names[imod])

        # Generate the model_list string
        model_list = '_model_list = ['
        for imod in range(len(mod_names)-1):
                model_list = model_list + ('"%s",'%mod_names[imod])

        model_list = model_list + ('"%s"'%mod_names[len(mod_names)-1])
        model_list = model_list + ']'
        wfile.write("%s\n"%model_list)

    print(' Modellist has been updated successfully.')
    return 


def getModelNames():
    """Returns the name of the available models.

    """

    mod_names = []
   
    # Get the name of all model files in the module directory
    import radmc3dPy
    mod_path = radmc3dPy.__file__.strip()[:-12]
    dum = sp.Popen(['ls -1 '+mod_path+'/models/*.py'], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE).communicate()[0].split()

    for i in range(len(dum)):
        sdum = str(dum[i].decode('utf-8')).split('/')
        modname = sdum[len(sdum)-1][:-3]

        if (modname!='template')&(modname!='__init__')&(modname!='_libfunc')&(modname!='_modellist'):
            mod_names.append(modname)

    return mod_names
    

def getModelDesc(model=None):
    """Returns a brief description of the selected model.
    
    Parameters
    ----------
    model   : str   
              Name of the model to get the description of

    Returns
    -------
    A string with the description of the model
    """

    if model is None:
        raise ValueError('Unknown model. No model name is given.')

    try:
        mdl = importlib.import_module(model)
    except ImportError:
        try:
            mdl = importlib.import_module('radmc3dPy.models.'+model)
        except ImportError:
            msg = model+'.py could not be imported. The model files should either be in the '\
                  + ' current working directory or in the radmc3d python module directory'
            raise ImportError(msg)

    if callable(getattr(mdl, 'getModelDesc')):
        return mdl.getModelDesc()
    else:
        raise RuntimeError(model+'.py does not contain a getModelDesc() function.')





