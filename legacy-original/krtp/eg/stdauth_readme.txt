10 September 2012

Matching species names to a list of accepted names requires some automated
matching of Authority.  This field involves punctuation and abbreviations of
people's names, and formats are not entirely standardized.  This issue became
especially obvious when merging Solanaceae Source and The Plant List
information.

My solution is to run any provided Authority string through a standardization
script.  There are two pieces to this script:

  * authority_alternates.dat
    This is a mini-text-database listing the various abbreviations employed for
    a given author.  If additional mis-matches are found, they can be added
    here.  I assembled the original contents by inspection of the Authority
    columns from the SS and PL data, by using google-refine, and with reference
    to an official(?) collection of names:
    http://www.ipni.org/ipni/authorsearchpage.do

  * std-auth.py
    This is the program that applies the standardization.  For input, it takes
    a file with one line per Authority string.  The output is the same, with
    each string standardized.

  * Typical usage:
    - Assemble a file called, e.g., "auth.in".  Perhaps make this by exporting
      a column from the table of names you have in hand to match.
    - Run the script: ./std-auth.py auth.in > auth.out  
    - Read "auth.out" into your table of names, replacing the original
      authorities column.
