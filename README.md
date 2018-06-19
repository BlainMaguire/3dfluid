# 3D Fluid

A 3D Real Time Fluid Solver based on Jos Stam's 2D fluid solver (stable Navier-Stokes solver).

[![Video Demonstration](http://img.youtube.com/vi/v4unokVKvrg/0.jpg)](https://www.youtube.com/watch?v=v4unokVKvrg "Video Demonstration")

Try it in your browser:

http://blainmaguire.com/projects/3dfluid/fluid.html

Reference: Jos Stam, "Real-Time Fluid Dynamics for Games". Proceedings of the Game Developer Conference, March 2003.

http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

### Controls:

* 'X' key - add source at center
* 'W' key - apply force x-axis
* 'D' key - apply force y-axis
* 'S' key - apply force z-axis
* 'C' key - clear simulation
* 'V' key - show/hide velocities
* 'A' Key - show/hide the XYZ axis
* 'H' key - show/hide this help message
* Left click  - pan from location
* Right click - rotate cube
* ESC - quit

### Potential Issues

**Nothing's Happening**

Be sure to press X a few times to add fluid. Then tap W, S, or D repeatedly. It could also be it's running way too slow on your computer.

**It runs really slow**

It takes exponentially more time to render larger the grid size gets. Try lowering the value of SIZE at the beginning of main.c. Increasing SIZE does make it look better but at the cost of performance.

**Compiling**

The code is all in c using GLUT, so it should be quite portable. I've included a simple makefile which I used.

The only thing I think might be problematic is I'm using GLUT to render bitmap characters for the help message via glutBitmapCharacter. This shouldn't really be a problem for newer implementations (openglut/freeglut).

If you are having other issues let me know.
