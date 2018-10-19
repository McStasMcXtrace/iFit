This class allows to communicate with a Matlab(R) instance, using a Pipe

Requirements:  python (2 or 3)
Optional:      ipython, ipython notebook, jupyter

Starting:
  First navigate to this location. You may as well just copy the matlab.py (and 
  optionally the PyMatlab.ipynb) file anywhere, and work from there.
  
  Then launch Python, e.g.:
    # ipython

  You can alternatively use a Python Notebook with:

    jupyter notebook PyMatlab.ipynb
  or 
    ipython notebook PyMatlab.ipynb

--------------------------------------------------------------------------------

The normal usage is, within Python:

>>> from matlab import Matlab
>>> m=Matlab()

By default the 'matlab' command is used, but you can specify an other program

>>> m=Matlab('/opt/MATLAB/R2010a/bin/matlab')
>>> m=Matlab('ifit')
    
Then you can get and set variables in the Matlab workspace
   
>>> m.set('a',1)
>>> m.get('a')
    
and directly evaluate any Matlab commands, including with control statements

>>> m.eval("a=1; disp('Now assigned variable a'); disp(a);")
>>> m.eval("if a, disp('a is true'); else disp('a is false'); end")

or even use multi-line strings and control statements (if, for, while, ...). 

By default, Python will wait for the commands to complete before returning 
to the prompt. The timeout is set with e.g.:

>>> m.timeout = 5

To force an asynchronous execution use:

>>> m.eval('commands',waitidle=False) 

especially for commands which are long to execute. In this case, use:

>>> m.waitidle()
>>> m.waitidle(timeout=10)

to wait for Matlab Idle state. You can stop this waiting with a keyboard
interrupt, leaving Matlab execution in asynchronous mode. 

Once you have finished your work with Matlab, you can close Matlab with:

>>> m.close()

which is also done automatically when exiting Python.
