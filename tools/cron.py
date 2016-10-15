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
Implements a Cron Scheduler in python.

Taken from stackoverflow.com:
http://stackoverflow.com/questions/373335/suggestions-for-a-cron-like-scheduler-in-python
"""

# ToDo: maybe add auto_start functionality of the CronDaemon (e.g. self.start within submit() method)?

import threading
import logging

import time
from datetime import datetime, timedelta

from tools.utility import Singleton, StoppableThread, timestamp

# Some utility classes / functions first
class AllMatch(set):
    """Universal set - match everything"""
    def __contains__(self, item): return True

allMatch = AllMatch()

def conv_to_set(obj):  # Allow single integer to be provided
    if isinstance(obj, (int,long)):
        return set([obj])  # Single item
    if not isinstance(obj, set):
        obj = set(obj)
    return obj

# The actual Event class
class CronEvent(object):
    def __init__(self, action, min=allMatch, hour=allMatch, 
                       day=allMatch, month=allMatch, dow=allMatch, 
                       args=(), kwargs={}):
        self.mins = conv_to_set(min)
        self.hours= conv_to_set(hour)
        self.days = conv_to_set(day)
        self.months = conv_to_set(month)
        self.dow = conv_to_set(dow)
        self.action = action
        self.args = args
        self.kwargs = kwargs

    def matchtime(self, t):
        """Return True if this event should trigger at the specified datetime"""
        return ((t.minute     in self.mins) and
                (t.hour       in self.hours) and
                (t.day        in self.days) and
                (t.month      in self.months) and
                (t.weekday()  in self.dow))

    def check(self, t):
        if self.matchtime(t):
            logging.getLogger().debug('Time match at cron event '+str(self)+' at time '+str(t)+'. Executing '+str(self.action)+'.')
            self.action(*self.args, **self.kwargs)

    def __repr__(self):
        return 'Cron Event on callable '+str(self.action)

class CronDaemon( ):
    __metaclass__ = Singleton
    
    def __init__(self, *events):
        self.events = list(events)
        self.lock = threading.Lock()
        self.thread = StoppableThread() # the thread the manager loop is running in

    def register(self, event):
        self.lock.acquire()
        try:
            self.events.append(event)
        finally:
            self.lock.release()

    def remove(self, event):
        self.lock.acquire()
        try:
            self.events.remove(event)
        finally:
            self.lock.release()

    def start(self):
        """Start the process loop in a thread."""
        if self.thread.is_alive():
            return
        self.thread = StoppableThread(target = self.run, name=self.__class__.__name__ + timestamp())
        self.thread.start()
 
    def stop(self, timeout=None):
        """Stop the process loop."""
        self.thread.stop(timeout=timeout)
 
    def run(self):
        logging.getLogger().info('Starting Cron Daemon.')
        t=datetime(*datetime.now().timetuple()[:5])
        while 1:
            
            self.lock.acquire()
            try:
                logging.getLogger().log(5,'Checking events at '+str(datetime.now())+' with '+str(t)+'.')
                for e in self.events:
                    e.check(t)
            finally:
                self.lock.release()
            
            t += timedelta(minutes=1)
            while datetime.now() < t:
                delta = (t - datetime.now()).total_seconds()
                self.thread.stop_request.wait(delta)
                if self.thread.stop_request.isSet():
                    return
        logging.getLogger().info('Shutting down Cron Daemon.')

if __name__ == '__main__':

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().info('Starting logger.')

    def act():
        print 'a'

    def bct():
        print 'b'

    CronDaemon().start()
    CronDaemon().register( CronEvent(act, min=range(60)) )
    CronDaemon().register( CronEvent(bct, min=range(0,60,2)) )
    