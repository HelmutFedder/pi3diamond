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
There are several distinct ways to go through different NVs and perform
certain measurement tasks. 

1. using the queue and 'SwitchTarget' and 'SaveJob' and 'SetJob'

  For each task, create a job and submit it to the queue.
  Provide a 'special' job for switching the NV. I.e., a queue might
  look like this: [ ODMR, Rabi, SwitchTarget, ODMR, Rabi, SwitchTarget, ...]

  pro: - very simple
       - a different set of Jobs can be submitted for individual NVs
       - every part of the 'code' is basically tested separately (uses only
         existing jobs) --> very low chance for errors
       - queue can be modified by user on run time, e.g., if an error in the tasks
         is discovered, it can be corrected
       - the submitted jobs can be run with lower priority than all the usual
         jobs, i.e., the queue can be kept during daily business and will
         automatically resume during any free time
         
  con: - no complicated decision making on how subsequent tasks are executed,
         e.g., no possibility to do first a coarse ESR, then decide in which range
         to do a finer ESR, etc. 
       - it is easy to forget save jobs. If everything goes well this is not a problem,
         because the jobs can be saved later at any time, but if there is a crash,
         unsaved jobs are lost

2. using an independent MissionControl job that is not managed by the JobManager

  Write a new job, that is not managed by the JobManager, i.e., that runs independently
  of the queue. This Job will submit jobs to the queue as needed.
  
  pro: - allows complex ways to submit jobs, e.g., depending on the result of previous
         measurement, with analysis performed in between, etc.

  con: - cannot be changed after started
       - control job will often be 'new code' and thus may have errors. It is
         difficult to test --> error prone
"""

from tools.emod import ManagedJob, JobManager

from tools.utility import timestamp

class SwitchTarget( ManagedJob ):
    
    def __init__(self, name=None):
        super(SwitchTarget, self).__init__()
        self.name=name
    
    def _run(self):
        
        try:
            self.state='run'
            if self.name is None:
                auto_focus.next_target()
            else:
                auto_focus.current_target=self.name            
            auto_focus._run()
        finally:
            self.state='idle'
            
    def __repr__(self):
        if self.name is None:
            return 'Switch to next Target'
        else:
            return 'Switch to '+str(self.name)

class SaveJob( ManagedJob ):
    
    def __init__(self, job, filename=None, timestamp=False):
        super(SaveJob, self).__init__()
        self.job = job
        self.filename = filename
        self.timestamp = timestamp
    
    def _run(self):
        
        try:
            self.state='run'
            job = self.job
            if self.filename is None:
                filename=str(job)
            else:
                filename=self.filename
            if self.timestamp:
                filename = timestamp() + '_' + filename
            job.save(filename)
        finally:
            self.state='idle'
            
    def __repr__(self):
        if self.filename is None:
            filename = str(self.job)
        else:
            filename = self.filename
        if self.timestamp:
            filename = '<timestamp>_' + filename
        return 'Save '+str(self.job)+' to '+filename

class SetJob( ManagedJob ):
    
    def __init__(self, job, d):
        super(SetJob, self).__init__()
        self.job = job
        self.d = d
    
    def _run(self):
        
        try:
            self.state='run'
            self.job.set_items(self.d)
        finally:
            self.state='idle'
            
    def __repr__(self):
        return 'Set items '+str(self.d)+' of job '+str(self.job)

for target in auto_focus.targets.keys():
    JobManager().submit(SwitchTarget(target))
    odmr = me.odmr.ODMR()
    odmr.stop_time=360
    JobManager().submit(odmr)
    JobManager().submit(SaveJob(odmr,target+'_odmr.pys',timestamp=True))
    