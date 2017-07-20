"""
A class to communicate with Matlab from a Python session.

Use matlab.Matlab? to get more help about how to use this module.
"""
from __future__ import print_function

__docformat__ = "restructuredtext en"
__version__   = '1.0'
__author__    = "E. Farhi <farhi@ill.fr>"

# all dependencies
import numpy
import scipy.io     # for MAT files I/O
import tempfile     # to create temporary .mat files
import atexit       # to close Matlab when exiting python
import subprocess   # for communication
import time         # for sleep at open
import sys
import os

# we use a Thread and a Queue to communicate asynchronously with matlab.
# indeed, the stdio may come any time. A separate thread then collects lines
# as they come, and make these available to the Matlab object when requested.
from threading  import Thread
try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty  # python 3.x
    
# a function to concatenate a queue from a stream, line by line
# this is used by a Thread as stdout reader
def enqueue_output(out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()
    
class MatlabError(Exception):
    """Raised when a Matlab evaluation results in an error inside Matlab."""
    pass

# ==============================================================================
class Matlab(object):
    """
    This class allows to communicate with a Matlab(R) instance, using a Pipe
    
    The normal usage is:
    
    >>> from matlab import Matlab
    >>> m=Matlab()
    >>> m.open()
    
    By default the 'matlab' command is used, but you can specify an other program
    
    >>> m.open('/opt/MATLAB/R2010a/bin/matlab')
    >>> m.open('ifit')
        
    Then you can get and set variables in the Matlab workspace
       
    >>> m.set('a',1)
    >>> m.get('a')
        
    and directly evaluate any Matlab commands, including with control statements
    
    >>> m.eval("a=1; disp('Now assigned variable a'); disp(a);")
    >>> m.eval("if a, disp('a is true'); else disp('a is false'); end")
    
    or even use multi-line strings and control statements (if, for, while, ...). 
    
    By default, Python will wait for the commands to complete before returning 
    to the prompt. 
    To force an asynchronous execution use:
    
    >>> m.eval('commands',waitidle=False) 
    
    especially for commands which are long to execute. In this case, use:
    
    >>> m.waitidle()
    
    to wait for Matlab Idle state. You can stop this waiting with a keyboard
    interrupt, leaving Matlab execution in asynchronous mode. 
    
    Once you have finished your work with Matlab, you can close Matlab with:
    
    >>> m.close()
    
    """
    

    def __init__(self):
        self.proc   = None  # the subprocess
        self.reader = None
        self.busy   = False # indicates when Matlab is ready/idle
        self.stdout = ''    # the stdout from Matlab
        self.executable = ''
        self.format = 'mat'
        
        # prompt: this is sent after every command to make sure we
        # detect the Idle state
        self.prompt = '>>'; 
        # error: displayed by Matlab when printing an error message
        self.error  = '???'

    def open(self, executable='matlab'):
        """
        Open a Matlab session
        input:
            executable: name of the Matlab executable (default is 'matlab')
        """
        
        self.executable = executable
        
        
        # we shall use stdout and stdin. Stderr is left untouched.
        # The os.setsid() is passed in the argument preexec_fn so
        # it's run after the fork() and before  exec() to run the shell.
        
        os.setpgrp() # create new process group, become its leader
        
        self.proc = subprocess.Popen([executable,'-nosplash','-nodesktop'], 
            preexec_fn=os.setsid,
            stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            
        # create an independent stdout reader sending data in a queue
        self.queue = Queue()
        self.thread= Thread(target=enqueue_output, 
                            args=(self.proc.stdout, self.queue))
        self.thread.daemon = True # thread dies with the program
        self.thread.start()
        print('Opened session', ' '.join([executable,'-nodesktop','-nosplash']), 'as PID', self.proc.pid)
        
        # make sure we close the pipe in all cases
        atexit.register(lambda handle=self: handle.close())
        
        # print start message reading stdout
        self.eval(waitidle=True) # check initial state
        self.busy = False
        
        return self
        # end open

    def close(self):
        """Close the Matlab session"""

        print('Closing session', self.executable, 'as PID', self.proc.pid)
        
        # force stop
        self.eval("exit",waitidle=False)  # from Matlab, gently
        time.sleep(.5)
        self.proc.stdin.close()           # emulates a Ctrl-D
        time.sleep(.5)
        self.proc.terminate()             # from outside, gently
        time.sleep(.5)
        self.proc.kill()                  # force
        self.eval() # force to poll the process and finalize close
        # end close
            
    # ==========================================================================
    def eval(self, expr=None, waitidle=True):
        """Evaluate an expression
        
        input:
            expr:     Matlab expression to evaluate
            waitidle: When True (default), execution is synchronous and Python
                      waits for the command to complete. When False, the Python
                      interpreter does not wait for Matlab and execution is
                      asynchronous. Use Matlab.waitidle() to wait.
        """
        
        # in this routine, at every eval call, we request to display the 
        # 'prompt' so that we can monitor when Matlab is ready.
        
        # check if the proc still exists
        if not self.proc:
            print('Matlab process is not started. Use Matlab.open()')
            return False
            
        self.proc.poll()
        
        # if no return code (still runs), we get its stdout
        if self.proc.returncode is not None:
            print('Matlab process is not active. Return code=', self.proc.returncode)
            return False
        
        if self.busy:
            print('eval: Matlab is busy. Use Matlab.waitidle() to wait for Idle state.')
            return
            
        if expr is not None and expr != '':
            # send the command
            self.proc.stdin.write(expr+"\n")            
        
        # request to display the prompt so that we can monitor the idle state
        self.proc.stdin.write("disp('"+self.prompt+"');\n")
        self.busy = True
        # check for idle state or timeout
        if waitidle:
            self.waitidle()     # flush until Idle
        else:
            busy = self.flush() # flush once
        # end eval
  
    def set(self, varname, value):
        """Set a variable name and assign it to the given value in Matlab"""

        # create a .mat file with the variable to send to matlab
        f = tempfile.NamedTemporaryFile()
        m = f.name+'.mat'
        scipy.io.savemat(m, { varname: value }, oned_as='row')
        
        # wait for file to be written
        time.sleep(.5)
        while not os.path.getsize(m):
            time.sleep(.5)
    
        # then call eval to read that variable within Matlab
        self.eval("load('"+m+"'); delete('"+m+"'); if isstruct("+varname+") && isfield("+varname+",'class'), try; "+varname+"=feval("+varname+".class, "+varname+"); end; end",
            waitidle=True)
        f.close()
        # end set
  
    def get(self, varname):
        """
        Get a variable from the Matlab workspace
        return:
            value:    requested variable
        """
        
        # check process, get its stdout and busy state
        # use eval to ask Matlab to save a variable into a .mat file
        f = tempfile.NamedTemporaryFile()
        m = f.name+'.mat'
            
        self.eval("tmp_class = class("+varname+"); if isobject("+varname+"), tmp_var=struct("+varname+"); tmp_var.class = tmp_class; else tmp_var="+varname+"; end; save('"+m+"', 'tmp_var','tmp_class','-v7'); clear tmp_var tmp_class;",
            waitidle=True)
            
        # wait for file to be written
        time.sleep(.5)
        while not os.path.getsize(m):
            time.sleep(.5)
    
        # then get that variable
        value = scipy.io.loadmat(m, squeeze_me=True, chars_as_strings=True)
            
        os.remove(m)
        f.close()
        
        # simplify object when loaded as multi-dimensional
        for key in value.iterkeys():
            if isinstance(value[key], numpy.ndarray):
                while value[key].shape and value[key].shape[-1] == 1:
                    value[key] = value[key][0]
        
        value = value['tmp_var']
        return value
        # end get
        
        
    # ==========================================================================
    def waitidle(self, timeout=60):
        """
        Wait for the Matlab interpreter to become idle
        
        input:
            timeout: a timeout in seconds for the wait. Default is 60.
        """
        
        start_time = time.time()
        
        # loop as long as in Busy state or timeout
        while self.flush():
            if timeout and time.time() - start_time > timeout:
                break
                
    def flush(self):
        """
        Flush the Matlab stdout stream and return the busy state.
        
        return:
            busy: False when Matlab is Idle, True when Busy.
        """
        
        # we flush all available stdout lines and wait for the prompt
        # the prompt is sent to Matlab stdin when calling eval.
        
        # read all available lines
        while True:
            try:
                line = self.queue.get_nowait()
            except Empty:
                # no more lines to read
                line = None
                break
            # got valid line.
            if line != '':
                # analyse the line, and search for the prompt (e.g. >>)
                rline = line.rsplit()
                if len(rline) > 0:
                    if rline[-1] == self.prompt:
                        self.busy = False
                    if rline[0] == self.error:
                        self.busy = False
                        raise MatlabError(line)
                    
                # store and display output
                line = line.replace(self.prompt,'')  # remove any prompt text for display
                self.stdout += line
                sys.stdout.write(line)
                sys.stdout.flush()

        return self.busy
        # end flush
        
# end Matlab class
