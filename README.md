compactnesslib
==============

compactnesslib is meant to be a performant library for calculating compactness
metrics of multipolygons, especially electoral districts.



Using it as a library
=====================

To use as a library, include the compactnesslib base directory in your project.

Include the library by adding:

    #include "compactnesslib/compactnesslib.hpp"

to your source file.

When you compile, compile all the `./*.cpp` and `./shapelib/*.c` files. The
files may be compiled all together without a need to respect the subdirectory
structure.



Available Scores/Metrics
========================

A list of the available compactness scores and how they are calculated is
available
[here](https://github.com/r-barnes/compactnesslib/blob/master/Scores.md).



Assumptions
========================

compactnesslib assumes that you have appropriately projected your data to a
Euclidean plane prior to performing any calculations on it.



Who uses compactnesslib?
========================

compactnesslib is a core component of the Python package
[python-mander](https://github.com/gerrymandr/python-mander)
and the R package
[mandeR](https://github.com/r-barnes/mandeR).



Credits
=======

This package was created for the
[Metric Geometry And Gerrymandering Group](https://sites.tufts.edu/gerrymandr/)
as part of a hack-a-thon on August 10-11, 2017.