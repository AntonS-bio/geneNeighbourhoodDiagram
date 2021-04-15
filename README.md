Requires: python and R
R needs to have gggenes package installed.

To run: python generateGeneNeighbourhoodDiagram.py AHM87981.1 3 ~/myGffDir/ MyGeneDiagram.jpg

AHM87981.1 - replace with your own genes ID, but should also work with any other type of ID. Script looks for this string in the gff line. 
3 - # of genes left and right of the above genes
~/myGffDir/ - where all gff files are
MyGeneDiagram.jpg - name of output file (optional)

Two parameters hardcoded in script for the moment (both at the top of generateGeneNeighbourhoodDiagram.py):
GffsToExclude - self explanatory. Values is list can be with and without ".gff" at the end.
minComboFrequency - any gene profile which occurs if fewer cases than this number will not be reported. By default 2. This is useful to remove occasional misassemblies if dataset is large.