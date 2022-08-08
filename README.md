**AbaqusFilFile-Translator** - a tool to convert .fil files produced by 
[Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) to
Ensight Gold .case files, which can be read by 
[EnSight](https://www.ansys.com/products/fluids/ansys-ensight) 
and [ParaView](https://www.paraview.org/).

This allows to perform postprocessing of Abaqus simulations without Abaqus/CAE.
In particular, results by simulations using UELs (user defined elements) can 
be postprocessed easily! No workaround using dummy elements in Abaqus/CAE is required 
anymore.

Requires Python 3.6+ (depends on dictionary with ordered storage) and Numpy.

Features
===========================
*  Custom conversion  of .fil files to Ensight .case files using export definition files.
*  Easy to use export definition file syntax.
*  Easy to extend! Written in Python, and currently unsupported .fil records can be added easily.
*  Conversion during Simulation! An ongoing Abaqus simulation is recognized, and the translator waits until new data is written by Abaqus.

Usage
===========================

    python filConverter.py file.fil exportDefinition.inp
    
Take a look at the example in the example directory.

Attention: Abaqus mangles the names of elsets and nodesets with part and assembly designations, and converts everything to uppercase. 
In order to identify your desired sets in the .fil file, just convert the .fil file with no exports defined (dry run) using the translator.
The translator will identify every existing set in the .fil file and print them in the console with the correct name!

Input keywords and options
===========================

general usage: 
--------------
    
    *keyword, option1=value1, option2=value2, ...
    **comments, blank lines are ignored

available keywords:
-------------------

    *UELSDVToQuadraturePoints    relate SDV data to quadrature points.

        destination                   string        new name of the result
        qpCount                       integer       define a periodical pattern for a repeatet extraction for results at
                                                    quadrature points
        qpDistance                    integer       define a periodical pattern: data distance between qps
        qpInitialOffset               integer       define a periodical pattern: initial constant offset before qp data
                                                    begins
        set                           string        Abaqus element set


    *computeAverageOverQuadraturePoints    perform a computation on an elemental result

        result                        string        Abaqus variable identifier
        set                           string        Abaqus element set


    *defineElementType    assign an ensight Shape to an Abaqus Element

        element                       string        Abaqus (User) Element
        shape                         string        Ensight Shape, can be any of: point bar2 bar3 tria3 tria6 quad4
                                                    quad8 tetra4 tetra10 pyramid5 pyramid13 penta6 penta15 hexa8 hexa20
                                                    nsided nfaced


    *ensightCaseOptions    modify Ensight export options

        discardTime                   string        discard Time values and replace by enumeration of time steps


    *ensightPerElementVariable    define an Ensight per element variable for export

        dimensions                    integer       (optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing
                                                    components will be zero filled
        exportName                    string        export name of the variable
        f(x)                          string        (optional), apply a mathematical/array expression on the result
                                                    array (per Element, slow!)
        location                      string        where is the result ? qps | computed
        result                        string        Abaqus variable identifier
        set                           string        Abaqus element set
        timeSet                       integer       (optional), define a timeset, for 'different' timelines
        values                        string        (optional), define a index/slice to extract a subarray from the
                                                    total result array (per Element)
        which                         string        which one? e.g. quadrature point numbers or "average" for average
                                                    computed results


    *ensightPerNodeVariable    define an Ensight per node variable for export

        dimensions                    integer       (optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing
                                                    components will be zero filled
        exportName                    string        export name of the variable
        f(x)                          string        (optional), apply a mathematical/array expression on the result
                                                    array (per Element, slow!)
        fillMissingValues             float         (optional), fill missing nodal values with a constant values,
                                                    requires specified dimensions (slow!)
        result                        string        Abaqus variable identifier
        set                           string        Abaqus element set
        timeSet                       integer       (optional), define a timeset, for 'different' timelines
        values                        string        (optional), define a index/slice to extract a subarray from the
                                                    total result array (per Element)


    *include    (optional) load extra .inp file (fragment), use relative path to current .inp

r       input                         string        filename
