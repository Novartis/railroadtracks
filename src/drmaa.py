# Copyright 2014 Novartis Institutes for Biomedical Research

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from railroadtracks.easy import (ParallelExecutionAbstract,
                                 _execute_subprocess_in_dir)
import itertools
import drmaa
import pyzmq
import threading


class ZMQServer(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        context = zmq.Context()
        frontend = context.socket(zmq.ROUTER)
        frontend.bind('')

class DRMAAExecution(ParallelExecutionAbstract):
    def __init__(self, processes):
        """ 
        :param processes: number of simultaneous Python processing (underlying tools may use threads or processes themselves)
        """
        self.__processes = processes

    def __get_ressources(self):
        pass
        
    def map(self, taskset, tempdir=None):
        assert isinstance(taskset, TaskSet)
        self.__get_ressources()
        # session
        session = drmaa.Session()
        session.initialize()
        # template
        jt = session.createJobTemplate() 
        jt.remoteCommand = sys.argv[1]
        jt.jobEnvironment = os.environ 
        jt.workingDirectory = sys.argv[2] 
        native_specs = ' '.join('-S /bin/bash',
                                '-l num_proc=%(processes)i' % {},
                                '-l interface=10gb',
                                '-pe orte 1',
                                ' -j y',
                                '-q default.q')
        jt.nativeSpecification = native_specs 
    
        command_args = "" 
        for i in range(2, len(sys.argv)): 
            command_args = command_args + str(sys.argv[i]) + ' ' 
            
        jt.args = [command_args] 
        jt.joinFiles=True 

        # working directories
        wds = (tempfile.mkdtemp(dir=tempdir) for x in taskset)
        # arguments for the tasks
        task_specs = ((x.task_id, x.unifex_cmd(), wd) for x, wd in zip(taskset, wds))
        tasklist = session.runBulkJobs(jt, _execute_subprocess_in_dir, task_specs)
        task.native_specification()
        
        # now we can start communicating with the remote workers
        context = zmq.Context()
        

        # Socket to receive messages on
        receiver = context.socket(zmq.PULL)
        receiver.connect("tcp://localhost:5557")

        # Socket to send messages to
        sender = context.socket(zmq.PUSH)
        sender.connect("tcp://localhost:5558")

        # Socket for control input
        controller = context.socket(zmq.SUB)
        controller.connect("tcp://localhost:5559")
        controller.setsockopt(zmq.SUBSCRIBE, b"")

        # Process messages from receiver and controller
        poller = zmq.Poller()
        poller.register(receiver, zmq.POLLIN)
        poller.register(controller, zmq.POLLIN)
        # Process messages from both sockets
        while True:
            socks = dict(poller.poll())

            if socks.get(receiver) == zmq.POLLIN:
                message = receiver.recv_string()

                # Process task
                workload = int(message)  # Workload in msecs

                # Do the work
                time.sleep(workload / 1000.0)

                # Send results to sink
                sender.send_string(message)

                # Simple progress indicator for the viewer
                sys.stdout.write(".")
                sys.stdout.flush()

            # Any waiting controller command acts as 'KILL'
            if socks.get(controller) == zmq.POLLIN:
                break


                project = next(iter(taskset)).project
                _set_task_status(project, res)
                for d in wds:
                    shutil.rmtree(d)


if __name__ == '__main__':
    pass
    # server
    
