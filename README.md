**AbaqusFilFile-Translator** - a tool to convert .fil files produced by 
[Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) to
Ensight Gold .case files, which can be read by 
[EnSight](https://www.ansys.com/products/fluids/ansys-ensight) 
and [ParaView](https://www.paraview.org/).

This allows to perform postprocessing of Abaqus simulations without Abaqus/CAE.
In particular, results by simulations using UELs (user defined elements) can 
be postprocessed easily! No workaround using dummy elements in Abaqus/CAE is required 
anymore.


Features
===========================
*  Custom conversion  of .fil files to Ensight .case files using export definition files
*  Easy to use export definition file syntax
*  Custom export of data to .csv files

Usage
===========================

    python filConverter.py <file>.fil exportDef.inp
    
Take a look at the example in the example directory.


Input keywords and options
===========================

general usage: 
--------------
    
    *keyword, option1=value1, option2=value2, ...
    dataline1_val1, dataline2_val1, dataline1_val3, ...
    dataline2 ...
    **comments, blank lines are ignored

available keywords:
-------------------
