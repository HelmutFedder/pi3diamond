# pi3diamond

pi3diamond is a toolkit for confocal scanning,
anti-bunching, FLIM, pulsed ODMR / NMR,
and more sophisticated quantum physics experiments,
typically performed with NV centers in diamond,
written in python using the enthought traits packages.

# Overview

pi3diamond is a set of python modules that consists of

* classes to control instruments (microwave sources, piezo scanners, etc.)
* widgets to run measurements (confocal scanning, ODMR, etc.)
* a scheduler and cron daemon to control the execution of measurements

To use the modules, you will typically write a little startup script,
that defines your hardware instruments and starts up GUI widgets
for the measurements you want to run.

# Quickstart

clone the git repository, copy 'example.py' to a new file, e.g. 'startup.py',
open it in an editor and adopt it to your needs.

open an ipython shell with a QT GUI event loop `ipython --gui=qt` and run
the startup script `run startup.py`

Note: you can try out most of pi3diamond without actual hardware connected
by using the mocking classes found in `hardware.mocking`.
