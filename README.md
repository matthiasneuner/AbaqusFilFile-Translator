Use abaqusFilConverter from console with:

python filConverter.py file.fil exportDef.inp

whith an exemplary exportDef.inp-file like e.g.:

*defineElementType, element=U3, shape=quad4
*ensightPerNodeVariable, set=mainPart, exportName=nodeDisplacements, source=U, dimensions=3,
0,1
*ensightPerElementVariable, set=ALL, exportName=ALPHAP_IP1, source=SDV
1
*ensightPerElementVariable, set=ALL, exportName=ALPHAD_IP1, source=SDV
3
*ensightPerElementVariable, set=ALL, exportName=STRESS_IP1, source=SDV
25,26,27,28,29,30
*ensightPerElementVariable, set=ALL, exportName=STRAIN_IP1, source=SDV
31,32,33,34,35,36