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

    *defineElementType    assign an ensight Shape to an Abaqus Element

        element                       string                        Abaqus (User) Element
        shape                         string                        Ensight Shape


    *ensightCaseOptions    modify Ensight export options

        discardTime                   string                        discard Time values and replace by enumeration of
                                                                    time steps


    *ensightPerElementVariable    define an Ensight per element variable for export

        dimensions                    integer                       (optional), 1/3/6/9 for
                                                                    scalar/vector/tensor/tensor9; missing components
                                                                    will be zero filled
        exportName                    string                        export name of the variable
        f(x)                          string                        (optional), apply a mathematical/array expression on
                                                                    the result array (per Element, slow!)
        integrationPointDataDistance  integer                       (optional), define a periodical pattern: initial
                                                                    constant offset )
        integrationPointDataOffset    integer                       (optional), define a periodical pattern: offset
                                                                    between extraction points
        nIntegrationPoints            integer                       (optional), define a periodical pattern for a
                                                                    repeatet extraction (e.g. for results @ GaussPts)
        set                           string                        Abaqus element set
        source                        string                        Abaqus variable identifier
        timeSet                       integer                       (optional), define a timeset, for 'different'
                                                                    timelines
        values                        string                        (optional), define a index/slice to extract a
                                                                    subarray from the total result array (per Element)


    *ensightPerNodeVariable    define an Ensight per node variable for export

        dimensions                    integer                       (optional), 1/3/6/9 for
                                                                    scalar/vector/tensor/tensor9; missing components
                                                                    will be zero filled
        exportName                    string                        export name of the variable
        f(x)                          string                        (optional), apply a mathematical/array expression on
                                                                    the result array (per Element, slow!)
        fillMissingValues             float                         (optional), fill missing nodal values with a
                                                                    constant values, requires specified dimensions
                                                                    (slow!)
        set                           string                        Abaqus element set
        source                        string                        Abaqus variable identifier
        timeSet                       integer                       (optional), define a timeset, for 'different'
                                                                    timelines
        values                        string                        (optional), define a index/slice to extract a
                                                                    subarray from the total result array (per Element)


    *include    (optional) load extra .inp file (fragment), use relative path to current .inp

        input                         string                        filename

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
