# BioDT

Do all your bioinformatics from R, using `data.table`s.

**WORK IS UNDERWAY converting our lab's raggedy but much beloved set of scripts into a more aerodynamic package suitable for public consumption.**

Represents common bioinformatics formats like fasta, bed, vcf, as `data.table`s with standardised columns, but which can be decorated or manipulated to the user's desire.

*Includes:*
- Read/write functions (file ==> data.table; data.table ==> file) and converters (data.table ==> data.table)
- Wrappers for common linux bioinformatics tools like blast and bedtools
- Interfaces e.g. to pluck out chunks of [bv]cf files as data.tables, using syntax like `vcf.interface["chromosome",start,end]`
- Various plotting and visualisation functions
- Various convenience functions for bioinformatics tasks that keep cropping up again and again (`alleleFrequency()`, `suggestKaspMarkers()`, `setRefAlleleMostCommon()`, ...)

BioDT will be implemented and maintained by [Tim Rabanus-Wallace and Kelly Rodgers](https://safes.unimelb.edu.au/research/cropgem-lab#people) from the [CropGEM lab](https://safes.unimelb.edu.au/research/cropgem-lab) at the University of Melbourne.
