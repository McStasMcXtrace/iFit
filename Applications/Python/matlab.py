# comments
"""
docstrings"""


__docformat__ = "restructuredtext en"
__version__   = '1.0'
__author__    = "E. Farhi <farhi@ill.fr>"

# all dependencies
import numpy
import hdf5storage  # to exchange temporary files with matlab
import scipy.io as sio  # fallback solution when hdf5storage fails
import tempfile     # to create temporary .mat files
import atexit       # to close Matlab when exiting python
import subprocess   # for communication
import time         # for sleep at open
import sys
import os
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
    >>> m.eval("if a, disp('a is true')\nelse disp('a is false'); end")
    
    or even use multi-line strings.
    
    Once you have finished your work with Matlab, you can save the session, e.g.
    
    >>> m.eval("save")
    
    and close the Matlab program.
    
    >>> m.close()
    
    To lad back the matlab session again, use:
    
    >>> m.eval("load")
    
    """
    

    def __init__(self):
        self.proc   = None  # the subprocess
        self.reader = None
        self.busy   = False # indicates when Matlab is ready/idle
        self.stdout = ''    # the stdout from Matlab
        self.executable = ''
        self.prompt = '>>'; # this is sent after every command to make sure we
                            # detect the Idle state

    def open(self, executable='matlab'):
        """
        Open a Matlab session
        input:
            executable: name of the Matlab executable (default is 'matlab')
        """
        
        self.executable = executable
        print 'Opening session ', ' '.join([executable,'-nodesktop','-nosplash'])
        
        # we shall use stdout and stdin. Stderr is left untouched.
        self.proc = subprocess.Popen([executable,'-nosplash','-nodesktop'], 
            shell=True,
            stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            
        # create an independent stdout reader sending data in a queue
        self.queue = Queue()
        self.thread= Thread(target=enqueue_output, 
                            args=(self.proc.stdout, self.queue))
        self.thread.daemon = True # thread dies with the program
        self.thread.start()
        
        # make sure we close the pipe in all cases
        atexit.register(lambda handle=self.proc: handle.kill())
        
        # print start message
        time.sleep(1)
        self.isbusy()  # check initial state

    def close(self):
        """Close the Matlab session"""
        
          
        # when Idle, we just do a terminate, else we kill.
        print 'Closing session ', self.executable
        
        # check process, get its stdout and busy state
        if self.isbusy():
            self.proc.kill()
        else:
            self.proc.terminate()
  
    def eval(self, expr=None, waitidle=False):
        """Evaluate an expression
        
        input:
            expr:     Matlab expression to evaluate
            waitidle: flag to wait for Matlab to become idle. Default is False.
            
        """
        
        # check process, get its stdout and busy state
        if self.isbusy():
            print 'eval: Matlab is busy.'
        elif expr is not None:
            self.proc.stdin.write(expr+"\n")
            if waitidle:
                self.waitidle()
  
    def set(self, varname, value):
        """Set a variable name and assign it to the given value in Matlab"""

        # create a scipy.io.savemat with the variable to send to matlab
        f = tempfile.NamedTemporaryFile()
        m = f.name+'.mat'
        
        # hdf5storage requires unicode dict members.
        # https://github.com/frejanordsiek/hdf5storage/issues/17
        if isinstance(value ,dict):
            value = dict([(k.decode(), v) for k, v in value.items()])
        
        try:
            hdf5storage.savemat(m, { varname: value }, 
                format='7.3', oned_as='column', store_python_metadata=True)
        except:
            hdf5storage.savemat(m, { varname: value }, 
                format='5', oned_as='column', store_python_metadata=True)
        
        # wait for file to be written
        time.sleep(.5)
        while not os.path.getsize(m):
            time.sleep(.5)
    
        # then call eval to read that variable
        self.eval("load('"+m+"'); delete('"+m+"');", waitidle=True)
        f.close()
  
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
        self.eval("tmp_class = class("+varname+"); if isobject("+varname+"), tmp_var=struct("+varname+"); tmp_var.class = tmp_class; else tmp_var="+varname+"; end; save('"+m+"', 'tmp_var','tmp_class', '-v7'); clear tmp_var tmp_class;",
            waitidle=True)
        
        # wait for file to be written
        time.sleep(.5)
        while not os.path.getsize(m):
            time.sleep(.5)
    
        # then get that variable
        try:
            value = hdf5storage.loadmat(m, variable_names='tmp_var')
        except IOError:
            value = sio.loadmat(m, variable_names='tmp_var',
                matlab_compatible=True, squeeze_me=True)
            value = value['tmp_var']
        os.remove(m)
        f.close()
        return value
        
    def waitidle(self, timeout=60):
        """
        Wait for the Matlab interpreter to become idle
        
        input:
            timeout: a timeout in seconds for the wait. Default is 60.
        """
        
        start_time = time.time()
        
        while self.isbusy():
            if timeout and time.time() - start_time > timeout:
                break
                
    def flush(self):
        """
        Flush the Matlab stdout stream and return the last line.
        """
        
        # read all available lines, and store the last one
        lastline = ''
        while True:
            try:
                line = self.queue.get_nowait()
            except Empty:
                # no more lines to read
                line = None
                break
            else: # got line.
                if line != '':
                    self.stdout += line
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    lastline = line
        return lastline

    def isbusy(self, waitidle=False):
        """
        Read the Matlab output and return the busy state.
        
        input:
            waitidle: flag to wait for Matlab to become idle. Default is False.
        """
        
        # a stdout asynchronous reader must be done inside a separate Thread
        # which only gets lines when they come. Then the main Thread (here)
        # will just poll if lines are available without blocking.
        #  see https://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python
        
        
        # check if the proc still exists
        self.proc.poll()
        
        # if no return code (still runs), we get its stdout
        if self.proc.returncode is not None:
            print 'Matlab process has died. Return code=', self.proc.returncode
            return False
            
        # send the prompt to check if Matlab responds
        self.proc.stdin.write("disp('"+self.prompt+"');\n")
        
        flag = True
        
        while flag:
            # read all available lines, and store the last one
            lastline = self.flush()
            
            # check the last line
            rlastline = lastline.rsplit()
            if len(rlastline) > 0 and lastline.rsplit()[-1] == self.prompt:
                busy     = False
                flag     = False
            else:
                busy = True
            self.busy = busy
            if not waitidle:
                break

        # store state and return it
        return self.busy
   
