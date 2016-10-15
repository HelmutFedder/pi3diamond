"""
    This file is part of pi3diamond, a toolkit for
    confocal scanning, anti-bunching, FLIM, pulsed ODMR / NMR,
    and more sophisticated quantum physics experiments,
    typically performed with NV centers in diamond,
    written in python using the enthought traits packages.

    pi3diamond is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pi3diamond is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diamond. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2009-2016 Helmut Fedder <helmut@fedder.net>
"""

"""
The execution model.
"""

import threading
import logging

from tools.utility import Singleton, StoppableThread, timestamp

from traits.api import HasTraits, Instance, Enum, Range, Button
from traitsui.api import View, Item, HGroup

# ToDo: maybe add auto_start functionality of the JobManager (e.g. self.start within submit() method)?

class Job( HasTraits ):

    """
    Defines a job.
    
    Methods:
    
        start():        starts the job
        stop(timeout):  stops the job
        _run():         actual function that is run in a thread
        
    Data:
    
        priority:   priority of the job (used by a job manager to schedule the job)
        state:      shows the current state of the job, 'idle', 'run' or 'wait'
    
      In the current execution model, a job should be re-startable.
    I.e., when a job is stopped before it is finished, upon next
    start, the work should be continued e.g. previously acquired
    data should be kept and accumulated.
    
      It is the user's task to ensure that previous data is
    handled correctly and to decide when a job should be continued
    and when it should be restarted as a new measurement. 

      A job can be in one of three states 'idle': doing nothing,
    'run': running, 'wait': waiting to be executed. The latter state
    is typically set by a Job manager to show that the job is
    scheduled for execution. The  
    """

    priority = Range(low=0, high=10, value=0, desc='priority of the job', label='priority', mode='text', auto_set=False, enter_set=True)

    state = Enum('idle', 'run', 'wait', 'done', 'error') # only for display. Shows the state of the job. 'idle': not submitted, 'run': running, 'wait':in queue

    thread = Instance( StoppableThread, factory=StoppableThread )

    def start(self):
        """Start the thread."""
        if self.thread.is_alive():
            return
        self.thread = StoppableThread(target = self._run, name=self.__class__.__name__ + timestamp())
        self.thread.start()
        
    def stop(self, timeout=None):
        """Stop the thread."""
        self.thread.stop(timeout=timeout)

    def _run(self):
        """Method that is run in a thread."""
        try:
            self.state='run'
            while(True):
                #logging.getLogger().debug("Yeah, still taking data like the LHC!")
                self.thread.stop_request.wait(1.0) # little trick to have a long (1 s) refresh interval but still react immediately to a stop request
                if self.thread.stop_request.isSet():
                    logging.getLogger().debug('Received stop signal. Returning from thread.')
                    break
            if True:
                self.state='idle'
            else:
                self.state='done'
        except:
            logging.getLogger().exception('Error in job.')
            self.state='error'
        finally:
            logging.getLogger().debug('Turning off all instruments.')            

class JobManager( ): # ToDo: In principle this need not be a singleton. Then there could be different job managers handling different sets of resources. However currently we need singleton since the JobManager is called explicitly on ManagedJob class.
    __metaclass__ = Singleton
    
    """Provides a queue for starting and stopping jobs according to their priority."""
        
    def __init__(self):
        self.thread = StoppableThread() # the thread the manager loop is running in
        self.lock = threading.Condition() # lock to control access to 'queue' and 'running'
        self.queue = []
        self.running = None
        self.refresh_interval = 0.1 # seconds
    
    def submit(self, job):

        """
        Submit a job.
        
        If there is no job running, the job is appended to the queue.

        If the job is the running job or the job is already in the queue, do nothing.
        
        If job.priority =< priority of the running job,
            the job is appended to the queue and the queue sorted according to priority.
        
        If job.priority > priority of the running job,
            the job is inserted at the first position of the queue, the running job is stopped
            and inserted again at the first position of the queue.
        """

        logging.debug('Attempt to submit job '+str(job))
        self.lock.acquire()
        
        running = self.running
        queue = self.queue

        if job is running or job in queue:
            logging.info('The job '+str(job)+' is already running or in the queue.')
            self.lock.release()
            return

        queue.append(job)
        queue.sort(cmp=lambda x,y: cmp(x.priority,y.priority), reverse=True) # ToDo: Job sorting not thoroughly tested
        job.state='wait'
                    
        logging.debug('Notifying process thread.')
        self.lock.notify()
            
        self.lock.release()
        logging.debug('Job '+str(job)+' submitted.')
 
    def remove(self, job):
        
        """
        Remove a job.
        
        If the job is running, stop it.
        
        If the job is in the queue, remove it.
        
        If the job is not found, this will result in an exception.
        """
 
        logging.debug('Attempt to remove job '+str(job))
        self.lock.acquire()

        try:
            if job is self.running:
                logging.debug('Job '+str(job)+' is running. Attempt stop.')
                job.stop()
                logging.debug('Job '+str(job)+' removed.')
            else:
                if not job in self.queue:
                    logging.debug('Job '+str(job)+' neither running nor in queue. Returning.')
                else:
                    logging.debug('Job '+str(job)+' is in queue. Attempt remove.')
                    self.queue.remove(job)
                    logging.debug('Job '+str(job)+' removed.')
                    job.state='idle' # ToDo: improve handling of state. Move handling to Job?
        finally:
            self.lock.release()
        
    def start(self):
        """Start the process loop in a thread."""
        if self.thread.is_alive():
            return
        logging.getLogger().info('Starting Job Manager.')
        self.thread = StoppableThread(target = self._process, name=self.__class__.__name__ + timestamp())
        self.thread.start()
    
    def stop(self, timeout=None):
        """Stop the process loop."""
        self.thread.stop_request.set()
        self.lock.acquire()
        self.lock.notify()
        self.lock.release()        
        self.thread.stop(timeout=timeout)
    
    def _process(self):
        
        """
        The process loop.
        
        Use .start() and .stop() methods to start and stop processing of the queue.
        """
        
        while True:
            
            self.thread.stop_request.wait(self.refresh_interval)
            if self.thread.stop_request.isSet():
                break
            
            # ToDo: jobs can be in queue before process loop is started
            # what happens when manager is stopped while jobs are running?
            
            self.lock.acquire()
            if self.running is None:
                if self.queue == []:
                    logging.debug('No job running. No job in queue. Waiting for notification.')
                    self.lock.wait()
                    logging.debug('Caught notification.')
                    if self.thread.stop_request.isSet():
                        self.lock.release()        
                        break
                logging.debug('Attempt to fetch first job in queue.')
                self.running = self.queue.pop(0)
                logging.debug('Found job '+str(self.running)+'. Starting.')
                self.running.start()
            elif not self.running.thread.is_alive():
                logging.debug('Job '+str(self.running)+' stopped.')
                self.running=None
                if self.queue != []:
                    logging.debug('Attempt to fetch first job in queue.')
                    self.running = self.queue.pop(0)
                    logging.debug('Found job '+str(self.running)+'. Starting.')
                    self.running.start()
            elif self.queue != [] and self.queue[0].priority > self.running.priority:
                logging.debug('Found job '+str(self.queue[0])+' in queue with higher priority than running job. Attempt to stop running job.')            
                self.running.stop()
                if self.running.state != 'done':
                    logging.debug('Reinserting job '+str(self.running)+' in queue.')
                    self.queue.insert(0,self.running)
                    self.queue.sort(cmp=lambda x,y: cmp(x.priority,y.priority), reverse=True) # ToDo: Job sorting not thoroughly tested
                    self.running.state='wait'
                self.running = self.queue.pop(0)
                logging.debug('Found job '+str(self.running)+'. Starting.')
                self.running.start()                
            self.lock.release()        
        
        
class ManagedJob( Job ):

    """
    Job with methods and buttons that submit the job to the JobManager.
    
    Methods:
    
        submit():     submit the job to the JobManager.
        remove():     remove the job from the JobManager.
        
    Data:
        
        state:        shows the current state of the job, 'idle', 'run', 'wait' or 'error'
        
    GUI:
    
        submit_button:    calls submit()
        remove_button:    calls remove()
        
    """

    submit_button = Button(label='submit', desc='Submit the measurement to the job manager.')
    remove_button = Button(label='remove', desc='Remove the measurement from the job manager. Stop it if necessary.')
    
    def submit(self):
        """Submit the job to the JobManager."""
        JobManager().submit(self) # we just submit again and again. The JobManager takes care that duplicate submissions are ignored 

    def remove(self):
        """Remove the job from the JobManager. Stop job if necessary."""
        JobManager().remove(self) # we just try to remove. The JobManager takes care that unknown jobs are ignored 

    def _submit_button_fired(self):
        """React to submit button. Submit the Job."""
        self.submit() 

    def _remove_button_fired(self):
        """React to remove button. Remove the Job."""
        self.remove()

    traits_view=View(HGroup(Item('submit_button', show_label=False),
                            Item('remove_button', show_label=False),
                            Item('priority'),
                            Item('state', style='readonly'),),
                     resizable=True)

class FreeJob( Job ):

    """
    Job with buttons that start the job without the JobManager.
    
    GUI:
    
        start_button:    calls start()
        stop_button:     calls stop()
        
    """

    start_button = Button(label='start', desc='Starts the measurement.')
    stop_button  = Button(label='stop',  desc='Stops the measurement.')
    
    def _start_button_fired(self):
        """React to submit button. Submit the Job."""
        self.start() 

    def _stop_button_fired(self):
        """React to remove button. Remove the Job."""
        self.stop()

    traits_view=View(HGroup(Item('start_button', show_label=False),
                            Item('stop_button', show_label=False),
                            Item('priority'),
                            Item('state', style='readonly'),),
                     resizable=True)


# ToDo: 1. integration of Job with typical experiments (HasTraits) and JobManager (class 'ManagedJob' ?)
#       2. fix JobManager
#       3. CRON Daemon
#       4. provide changing of priority after submission

# ToDo: there is a fundamental problem with the execution model: a user could start to play with priorities of his measurements
#       and eventually give them priorities that are higher than the tracker priority --> this is a bit difficult to see

# ToDo: migrate Job to a subclass of Thread?



if __name__ == '__main__':
    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    JobManager().start()
    
    import time
    time.sleep(0.1)
    
    jobs = [ Job() for i in range(5)]
        
    jobs[0].priority = 0
    jobs[1].priority = 1
    jobs[2].priority = 0
    jobs[3].priority = 3
    jobs[4].priority = 0
    
    [JobManager().submit(job) for job in jobs ]

    time.sleep(0.1)
    q = JobManager().queue

    print [job.priority for job in q]
    print [q.index(job) if job in q else None for job in jobs]
    
    #import time
    #time.sleep(0.1)

    #j3 = Job()
    #j3.priority = 1
    #JobManager().submit(j3)
    
    