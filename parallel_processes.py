import subprocess
import time
import locale
import sys

os_encoding = sys.getfilesystemencoding()


def arg_list_to_str(args):
    return " ".join(args) # could be more sophisticated

class SimpleSubProcess:
    def __init__(self, args, number):
        self.process = subprocess.Popen(args, 
             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.returncode = None
        self.args = arg_list_to_str(args)
        assert number >= 0
        self.number = number

    def poll(self):
        if self.returncode is None:
             self.returncode = self.process.poll()
        return self.returncode    

    def wait(self):
        if self.returncode is None:
             self.returncode = self.process.wait()
        return self.returncode

    def kill(self):
         if self.returncode is None and self.process.poll() is None:
             self.process.kill()
         self.returncode = -1
         

    def display_output(self):
        self.wait()
        print(self.args) 
        output = self.process.stdout.read()
        try: 
            if len(output):
                print(output.decode(os_encoding, errors='backslashreplace'))
        except:
            print(output)


class SimpleQutputQueue:
      def __init__(self):
          self.display_next = 0
          self.processes = {}
      
      def display_finished(self):
          keys = sorted(self.processes.keys())
          for key in keys:
              if key <= self.display_next:
                  self.processes[key].display_output()
                  del self.processes[key]
                  self.display_next += 1

      def display_all(self):
          keylist = list(self.processes.keys()) + [self.display_next]
          self.display_next  = max(keylist)
          self.display_finished()
          self.reset()

      def enter_process(self, process):
          assert isinstance(process, SimpleSubProcess)
          num = process.number
          if num in self.processes:
              self.processes[num].display_output()
          self.processes[num] = process
          self.display_finished() 

      def reset(self):  
          self.display_next = 0
          self.processes.clear()



class SimpleProcessWorkers:
    def __init__(self, nprocesses = 1):
        self.nprocesses = max(nprocesses, 1)
        if self.nprocesses > 1:
            self.processes = [None] * self.nprocesses
            self.output_queue = SimpleQutputQueue()
            self.index = 0

    def _wait_input(self):
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
                process.display()
                self.processes[self.index] = None
                self.kill_processes()
                return -1

    def _wait_ready(self):
        for i in range(self.nprocesses):
            process = self.processes[self.index]
            if process is not None:
                status = process.wait()
                if status == 0:
                    self.output_queue.enter_process(process)
                    self.processes[self.index] = None
                else:
                    process.display()
                    self.processes[self.index] = None
                    self.kill_processes()
                    return -1
            self.index = (self.index + 1) % self.nprocesses      
        self.output_queue.display_all()
        self.output_queue.reset()
        self.index = 0
        return 0
           
                          
    def kill_processes(self):
        for i, process in enumerate(self.processes):
            if process is not None:
                 process.kill()  
            self.processes[i] = None
        self.output_queue.reset()
        self.index = 0

    def run_multiprocess(self, arglist):
        index = 0
        for i, args in enumerate(arglist):
            index = self._wait_input()
            if index < 0:
                break
            self.processes[index] = SimpleSubProcess(args, i)
        if index >= 0:
            index =  self._wait_ready()
        if index < 0:
            raise ValueError("Subprocess has failed")
           
    def run(self, arglist):
        if self.nprocesses > 1:
            self.run_multiprocess(arglist)
        else:
            for args in arglist:
                print(arg_list_to_str(args))
                subprocess.check_call(args) 

          
