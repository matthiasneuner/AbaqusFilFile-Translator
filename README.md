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

<p align="center">
  <img width="512" height="512" src="./share/hacc_disp_boomerang_10c.gif">
</p>


Features
===========================
*  Custom conversion  of .fil files to Ensight .case files using export definition files.
*  Easy to use export definition file syntax.
*  Easy to extend! Written in Python, and currently unsupported .fil records can be added easily.
*  Conversion during Simulation! An ongoing Abaqus simulation is recognized, and the translator waits until new data is written by Abaqus.

Usage
===========================

    python filconverter.py file.fil exportDefinition.inp
    
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

    *UELSDVToQuadraturePoints    relate SDV data to quadrature points

        destination                   string        new name of the result
        qpCount                       integer       define a periodical pattern for a repeated extraction for results at
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


    *ensightPerElementVariableJob    define an Ensight per element variable for export

        dimensions                    integer       (optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing
                                                    components will be zero filled
        name                          string        export name of the variable
        timeSet                       integer       (optional), define a timeset, for 'different' timelines


    *ensightPerElementVariableJobEntry    define an Ensight per element variable entry for an element set

        f(x)                          string        (optional), apply a mathematical/array expression on the result
                                                    array (per Element, slow!)
        job                           string        export name of the variable
        location                      string        where is the result ? qps | computed
        result                        string        Abaqus variable identifier
        set                           string        Abaqus element set
        values                        string        (optional), define a index/slice to extract a subarray from the
                                                    total result array (per Element)
        which                         string        which one? e.g. quadrature point numbers or "average" for average
                                                    computed results


    *ensightPerNodeVariableJob    define an Ensight per node variable for export

        dimensions                    integer       (optional), 1/3/6/9 for scalar/vector/tensor/tensor asym; missing
                                                    components will be zero filled
        name                          string        export name of the variable
        timeSet                       integer       (optional), define a timeset, for 'different' timelines


    *ensightPerNodeVariableJobEntry    define an Ensight per node variable for an element set

        f(x)                          string        (optional), apply a mathematical/array expression on the result
                                                    array (per Element, slow!)
        fillMissingValues             float         (optional), fill missing nodal values with a constant values,
                                                    requires specified dimensions (slow!)
        job                           string        The associated export job
        result                        string        Abaqus variable identifier
        set                           string        Abaqus setname
        setType                       string        elSet or nSet, default=elSet
        values                        string        (optional), define a index/slice to extract a subarray from the
                                                    total result array (per Element)


    *ignoreLastNodesForElementType    Ignore trailing nodes to be ignored (e.g, make a hexa27 to a hex20 with number=7)

        element                       string        Abaqus (User) Element
        number                        integer       The number of nodes to be ignored


    *include    (optional) load extra .inp file (fragment), use relative path to current .inp

        input                         string        filename


    *substituteElSet    define an substitution for an element set in the .fil file. This is useful if you want to
                        replace an element set with another one.For instance Abaqus/Explicit is know to write faulty
                        element sets to the *.inp file if multiple cpu cores are used in combination with VUEL/VUMAT.

        data                          string        Abaqus like element set definition lines, i.e., the list of element
                                                    labels.
        elSet                         string        The name of the set to be substituted
