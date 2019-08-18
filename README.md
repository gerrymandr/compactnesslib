compactnesslib
==============

compactnesslib is meant to be a performant library for calculating compactness
metrics of multipolygons, especially electoral districts.



Using it as a library
=====================

To use as a library, include compactnesslib in your project's `CMakeLists.txt` file:

    add_subdirectory(compactnesslib)

Then include the library in your code using:

    #include <compactnesslib/compactnesslib.hpp>



Running Tests
=============

To run tests build using:

    mkdir build
    cd build
    cmake ..
    make unittest
    



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