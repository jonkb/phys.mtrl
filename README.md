# Phys.Mtrl: A Mechanics of Materials Workbench
phys.mtrl is GUI Workbench for performing Mechanics of Materials calculations.

## Dependencies
- Python 3
- tkinter
- numpy
- sympy
- matplotlib

## Implemented Features
Here is a summary of some of the currently implemented features.

### Interface
- Options
	- Set the scale of the lab
	- Toggle weightless members
- Simple interface for placing members of several cross sections and materials and placing supports and loads on those members
- File Save and Open: human-readable xml files are used to store the current state of the lab to be loaded later

### Evaluations
- Mass Properties: Report mass, volume, and several other values related to the member and its material
- Axial Stress and Strain: Reports the stress and strain caused only by axial forces (compression and tension) on the member
- Shear and Moment: Plots Shear and Moment along the length of a beam
- Axial and Shear Stress: Plots Axial and Shear stresses at every point in the profile of the member, based on axial forces combined with bending. Also calculates the principle stresses and identifies max and min values
- Euler Buckling: Reports if the column is stable according to the Euler Buckling theory

## Possible Future Improvements
- Make it so you can click anywhere on the member to select it, not just near the axis. This matters for exceptionally wide members.
- Finish the Euler Buckling evaluation. Currently it only works for fixed-free columns with one load exactly at the end.
- Rework the system to allow members to be placed at angles instead of just horizontal and vertical.
- Allow members to be moved after being placed.
- Allow forces and supprts to be deleted after being placed.
- Add joints to connect members together.
- Use equations for elongation to deal with statically indeterminate members.
- Add torque.
- Add an option for snap to grid when placing everything.
- Add more options for cross section and material.
	- Upload some database of materials and material properties
- Fix the mohr\_transform function so it works fine with symbolic functions, then use that to find the max and min symbolically in the sig\_tau report.
- Make it so that the max2d\_grad function works better with piecewise.
- Find a better way to represent the principal stresses when they're only known along the neutral axis, like for a circular cross section (sig\_tau report).
- Nicer support appearance, and supports that can be placed anywhere along the axis.
- Add keyboard shortcuts

## LICENSE
MIT License. See LICENSE file.


