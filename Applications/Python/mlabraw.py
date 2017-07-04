#!/usr/bin/env python

""" A quick and extremely dirty hack to wrap matlabpipe/matlabcom as if they
were mlabraw.

Author: Dani Valevski <daniva@gmail.com>
License: MIT
"""
import sys

is_win = sys.platform.startswith('win')
if is_win:
  from matlabcom import MatlabCom as MatlabConnection
  from matlabcom import MatlabError as error
else:
  from matlabpipe import MatlabPipe as MatlabConnection
  from matlabpipe import MatlabError as error

try:
  import settings
except:
  class settings:
    MATLAB_PATH = 'guess'

def open(matlab_path=None):
  if matlab_path is None:
    # check if we can find Matlab
    if os.getenv('MLABRAW_CMD_STR'):
      matlab_path = os.getenv('MLABRAW_CMD_STR')
    else:
      matlab_path = settings.MATLAB_PATH
    if matlab_path != 'guess' and os.path.isfile(matlab_path + '/bin/matlab'):
      matlab_path = matlab_path + '/bin/matlab'
  if is_win:
    ret = MatlabConnection()
    ret.open()
  else:
    ret = MatlabConnection(matlab_process_path=matlab_path)
    ret.open()
    #except:
    #  print 'Could not open matlab, is it in %s?' % matlab_path
    # ret = None
  return ret
  
def close(matlab):
  matlab.close()

def eval(matlab, exp, log=False):
  if log or is_win:
    matlab.eval(exp)
  else:
    matlab.eval(exp, print_expression=False, on_new_output=sys.stdout.write)
  return ''

def get(matlab, var_name):
  return matlab.get(var_name)

def put(matlab, var_name, val):
  matlab.put({var_name : val})
