Version: $Revision$
$

These scripts require the iFunc class to be used. They create Models (for e.g. fitting)

Location:
--------------------------------------------------------------------------------
.             File in this directory, when called, directly return an iFunc Model.
Specialized   File in this directory, when called, directly return an iFunc Model
                for specialized use, e.g. neutron scattering.
Factory       File in this directory, when called, require additional parameters
                in order to create a Model. The ceation process may be long.
                Factory Models without input argument should display a Dialogue
                to enter creation parameters. 
                  Input arguments:
                  None or empty:  display a Dialogue
                  'identify':     return a default iFunc object to identify the Model
                  'defaults':     return a default iFunc object
                  other:          proceed with Model creation with other arguments



