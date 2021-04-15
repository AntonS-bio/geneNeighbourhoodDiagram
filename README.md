Requires: python and R<br />
R needs to have gggenes package installed.<br />

To run: python generateGeneNeighbourhoodDiagram.py AHM87981.1 3 ~/myGffDir/ MyGeneDiagram.jpg<br />
<br />
AHM87981.1 - replace with your own genes ID, but should also work with any other type of ID. Script looks for this string in the gff line.<br />
3 - # of genes left and right of the above genes<br />
~/myGffDir/ - where all gff files are<br />
MyGeneDiagram.jpg - name of output file (optional)<br />
 

Two parameters hardcoded in script for the moment (both at the top of generateGeneNeighbourhoodDiagram.py):<br />
<br />
GffsToExclude - self explanatory. Values is list can be with and without ".gff" at the end.<br />
minComboFrequency - any gene profile which occurs if fewer cases than this number will not be reported. By default 2. This is useful to remove occasional misassemblies if dataset is large.<br />
