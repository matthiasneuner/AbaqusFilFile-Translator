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

Usage
===========================

    python filConverter.py <file>.fil exportDef.inp
    
Take a look at the example in the example directory.


Input keywords and options
===========================

general usage: 
--------------
    
    *keyword, option1=value1, option2=value2, ...
    **comments, blank lines are ignored

available keywords:
-------------------

element type can be any of:
point g_point
bar2 g_bar2
bar3 g_bar3
tria3 g_tria3
tria6 g_tria6
quad4 g_quad4
quad8 g_quad8
tetra4 g_tetra4
tetra10 g_tetra10
pyramid5 g_pyramid5
pyramid13 g_pyramid13
penta6 g_penta6
penta15 g_penta15
hexa8 g_hexa8
hexa20 g_hexa20
nsided g_nsided
nfaced g_nfaced
