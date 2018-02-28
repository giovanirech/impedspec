# impedspec
Python Graphic User Interface for impedance spectroscopy analysis
-----------------------------------------------------------------------
|                              ImpedSpec                              |
|                                                                     |
|	 A Graphical User Interface for Impedance Spectroscopy using        |
|                             Keysight 4294A                          |
|																	                                    |
|                      (Version 3.1, October 2016)                    |
|																	                                    |
|          Copyright (C) 2018 G. L. Rech  and C. A. Perottoni         |
|                                                                     |
-----------------------------------------------------------------------
-----------------------------------------------------------------------

ImpedSpec is a graphical user interface to control of instrument, 
acquisition and processing of impedance spectra from dielectric materials 
and electric circuits.

This program is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
                                                                     
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
                                                                     
You should have received a copy of the GNU General Public License 
along with this library.  If not, see <http://www.gnu.org/licenses/>.

Publications reporting results obtained using ImpedSpec should make
appropriate citation to the software and papers that describe it.


-----------------------------------------------------------------------
								ImpedSpec
-----------------------------------------------------------------------

In this directory you will find the following files: 

script.py ---> The main python script. To run the software, you should run
				       run this file using: python script.py

interdace-ui.py ---> This is the user interface python file. This file is 
					     automatically generated from the .ui file using pyuic4.

InterfaceImpedSpec.ui ---> Qt file containing the user interface as developed
						   using QtDesigner. 

sample*.txt    ---> Example files from analysis performed using ImpedSpec. This
					     files can be visualized in the software Through File>Import Data.

LICENSE         ---> Copy of the GNU General Public License, Version 3

-----------------------------------------------------------------------
                            DEPENDENCIES
-----------------------------------------------------------------------
Notice that you should have the drivers of the equipment correctly installed
on your machine for the software to recognize it. Check the equipment support
page at Keysight website for the latest driver available.

Some python packages are needed (in addition to Numpy and matplotlib). You 
can get them by running

	pip install astropy docutils pyvisa pyvisa-sim

You should have Qt installed in your machine. Also, PyQt is needed.
You can install it through

	sudo dnf install pyqt4-devel
	
If you are having any trouble related to the matplotlib backend, 
try to run

	sudo dnf install python-matplotlib-qt
	
	
-----------------------------------------------------------------------
						RUNNING THE SOFTWARE
-----------------------------------------------------------------------
Once the dependencies are satisfied, the GUI can be excecuted using the command
prompt in the program folder :

	python script.py


-----------------------------------------------------------------------
						MODIFYING THE SOFTWARE
----------------------------------------------------------------------
If you wish to modify the software, you can do it by modyfing the script.py file
and by modifying the user interface. The GUI file (InterfaceImpedSpec.ui) can be
eddited using QtDesigner. In fedora you can get it by running

	sudo dnf install qt-devel qt-creator
	
After any modification in the GUI, you should create a new python interface file.
You can do it with pyuic4 (included with pyqt)

	pyuic4 InterfaceImpedSpec.ui -o interface_ui.py
	
	
