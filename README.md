# minfi-modified-functions

As described in poster 436F, ASHG 2015 annual meeting,
"Adjusting Infinium methylation profiles to suppress signals from varying cell proportion". 
The main content is a file containing a modified version of estimateCellCounts, and additional minfi functions called by 
this function. A second file is example code processing raw Illumina Infinium methylation data, using insertSource() 
to insert the modified functions into the minfi package, calling the modified estimateCellCounts, and processing the output 
from the modified function.  This estimateCellCounts output is a list with additional components than that returned by 
the standard estimateCellCounts, in particular an Mset adjusted for cell composition.   
