"""Support execution of subprogesses in parallel

Assume that we want to execute two simple subprocesses
with module subprocess, e.g.:

import subprocess
subprocess.check_call(['cc', 'hello1.c'])
subprocess.check_call(['cc', 'hello2.c'])  

This module supports parallel execution of these two processes
as follows:

from parallel_processes import SimpleProcessWorkers
args_list = [['cc', 'hello1.c'], ['cc', 'hello2.c']]
workers = SimpleProcessWorkers(2) # up to 2 parallel processes
workers.run(args_list)

"""


import subprocess
import time
import locale
import sys
from subprocess import CalledProcessError

os_encoding = sys.getfilesystemencoding()


def arg_list_to_str(args):
    """Map list of arguments of a subbprcess to a string"""
    return " ".join(args) # could be more sophisticated

class SimpleSubProcess:
    """Models a simple subprocess"""
    PIPESIZE = 0x80000
    def __init__(self, args, number):
        """Start a subprocess with argument list ``args``.

        Subprocesses should be started with a ``number``.
        Numbers should be consecutive, starting with 0.
        """
        self.returncode = None
        self.args = arg_list_to_str(args)
        try: 
            self.process = subprocess.Popen(args, 
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                 pipesize=self.PIPESIZE)
        except:
            # Some installations don't like the 'pipesize' keyword arg
            try:
                self.process = subprocess.Popen(args, 
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            except:
                sys.stdout.flush()
                sys.stderr.flush()
                ERR= """
Error: Could not launch the following subprocess!
 %s
"""
                print(ERR % self.args)
                raise
        assert number >= 0
        self.number = number

    def poll(self):
        """Poll a subprocess.

        The function returns ``None`` is the subprocess is still 
        runnning. If the subprocess has terminated the method returns
        returns its exit status as in method ``wait``.
        """
        if self.returncode is None:
             self.returncode = self.process.poll()
        return self.returncode    

    def wait(self):
        """Wait until subprocess has finished and return exit status

        Usually, the exit status is 0 if the process has terminated
        successfully, and a nonzero integer oherwise.
        """
        if self.returncode is None:
             self.returncode = self.process.wait()
        return self.returncode

    def kill(self):
        """Kill subprocess if it is still running"""
        if self.returncode is None and self.process.poll() is None:
            self.process.kill()
        self.returncode = -1
         

    def display_output(self):
        """Display output of subprocess to stdout.

        This function tries to display the output in the same way
        as if it had been started with the ``subprocess`` package.
        """
        self.wait()
        sys.stdout.flush()
        sys.stderr.flush()
        if self.returncode != 0:
             print("\nError: The following subprocess has failed!")
        print(self.args) 
        output = self.process.stdout.read()
        try: 
            if len(output):
                print(output.decode(os_encoding, errors='backslashreplace'))
        except:
            print(output)
        sys.stdout.flush()
        sys.stderr.flush()


class SimpleQutputQueue:
      """Queue of terminated proceseses waiting to diplay their output

      Here we try to display the output of all successfully terminated
      processs in the order in which they have been launched.

      Here process is is an instance of class ``SimpleSubProcess``;
      the order of the processes is given by attribute ``number``.     
      """
      def __init__(self):
          self.next = 0
          self.processes = {}
      
      def display_finished(self):
          """Dislay output of all processes with number <= self.next"""
          keys = sorted(self.processes.keys())
          for key in keys:
              if key <= self.next:
                  self.processes[key].display_output()
                  del self.processes[key]
                  self.next += 1

      def display_all(self):
          """Display output of all processes in the queue"""
          keylist = list(self.processes.keys()) + [self.next]
          self.next  = max(keylist)
          self.display_finished()
          self.reset()

      def enter_process(self, process):
          """Enter a terminated process into the queue

          ``process`` must be an instance of class ``SimpleSubProcess``.
          The also function displays all output in sequential order
          up to the point where a process has not yet terminated. 
          """
          assert isinstance(process, SimpleSubProcess)
          num = process.number
          if num in self.processes:
              self.processes[num].display_output()
          self.processes[num] = process
          self.display_finished() 

      def reset(self):
          """Reset the queue so that it will be empty"""  
          self.next = 0
          self.processes.clear()



class SimpleProcessWorkers:
    """Model a system of ``nprocesses`` workers for subprocesses

    An instance of this class can up to run ``nprocesses`` 
    subprocesses in parallel

    Here we support a rather simple model of subprocesses. Each
    subprocess returns with an exit status, which is 0 in case 
    of success and a non-zero integer in case of failure. 

    An arbitrary number of subprocesses can be launched with 
    method ``run``. 
    """
    def __init__(self, nprocesses = 1):
        self.nprocesses = max(nprocesses, 1)
        if self.nprocesses > 1:
            self.processes = [None] * self.nprocesses
            self.output_queue = SimpleQutputQueue()
            self.index = 0

    def _wait_input(self):
        """Wait until less then ``nprocesses`` are running

        The function returns the number ``i`` of a free slot such
        that the next process can be launched and memorized in 
        self.processes[i]. It returns -1 if a process has failed.
        Any terminated processes will be entered into the queue
        self.output_queue; and 
        """
        n_polled = 0
        while True:
            process = self.processes[self.index]
            if process is None:
                return self.index 
            status =  process.poll()
            if status is None:
                n_polled += 1
                if n_polled >= self.nprocesses:
                    time.sleep(0.05)
                    n_polled = 0
                self.index = (self.index + 1) % self.nprocesses      
            elif status == 0:
                self.output_queue.enter_process(process)
                self.processes[self.index] = None
                return self.index
            else:
                process.display_output()
                self.processes[self.index] = None
                self.kill_processes()
                return -1

    def _wait_ready(self):
        """Wait until all processes have terminated.

        The function returns 0 if all processes have terminated
        sucessfully. It stops and kills all other processes and
        returns -1 if any process has failed. 
        """
        for i in range(self.nprocesses):
            process = self.processes[self.index]
            if process is not None:
                status = process.wait()
                if status == 0:
                    self.output_queue.enter_process(process)
                    self.processes[self.index] = None
                else:
                    process.display_output()
                    self.processes[self.index] = None
                    self.kill_processes()
                    return -1
            self.index = (self.index + 1) % self.nprocesses      
        self.output_queue.display_all()
        self.output_queue.reset()
        self.index = 0
        return 0
           
                          
    def kill_processes(self):
        """Kill all subprocesses yet running"""
        for i, process in enumerate(self.processes):
            if process is not None:
                 try:
                     process.kill()  
                 except:
                     pass
            self.processes[i] = None
        self.output_queue.reset()
        self.index = 0

    def run_multiprocess(self, arglist):
        """Run subprocesses in parallel

        The function waits until all processes have terminated.
        Not for public use; use method ``run`` instead!
        """
        index = 0
        try:
            for i, args in enumerate(arglist):
                index = self._wait_input()
                if index < 0:
                    break
                #print("Launching", " ".join(args))
                self.processes[index] = SimpleSubProcess(args, i)
            if index >= 0:
                index =  self._wait_ready()
        except:
            self.kill_processes()
            raise
        if index < 0:      
            raise ValueError("Subprocess has failed")
           
    def run(self, arglist):
        """Run subprocesses in parallel

        Each subprocess is controlled ba a list of arguments as
        in the standard module ``subpocess``. Parameter ``arglist``
        is a list of such lists. For each entry of ``arglist`` a
        subprocess is launched with arguments given by that entry.

        The function waits until all subprocesses have terminated 
        successfuly, or one of the processes fails.

        The function raises ValueError or CalledProcessError
        if a process fails.
        """
        if self.nprocesses > 1:
            self.run_multiprocess(arglist)
        else:
            for args in arglist:
                print(arg_list_to_str(args))
                subprocess.check_call(args) 

          
