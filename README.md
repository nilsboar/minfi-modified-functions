# minfi-modified-functions

As described in poster 436F, ASHG 2015 annual meeting,
"Adjusting Infinium methylation profiles to suppress signals from varying cell proportion". 


The main content is a file containing a modified version of estimateCellCounts, and additional minfi functions called by 
this function. A second file is example code processing raw Illumina Infinium methylation data, using insertSource() 
to insert the modified functions into the minfi package, calling the modified estimateCellCounts, and processing the output 
from the modified function.  This estimateCellCounts output is a list with additional components than that returned by 
the standard estimateCellCounts, in particular an Mset adjusted for cell composition.   

The non-standard aspect of the methylation signal processing is that the users 450K data is corrected and normalized 
together with the FlowSorted.Blood.450k library data. This requires editing the pdata for the two sets to make them 
consistent so that the sets can be combined.
