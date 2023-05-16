## R CMD check results

0 errors | 0 warnings | 3 notes

* Checking installed package size:
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      doc   4.6Mb

* checking DESCRIPTION meta-information ... NOTE
  Maintainer field differs from that derived from Authors@R
    Maintainer: ‘Elizabeth Holmes - NOAA Federal <eli.holmes@noaa.gov>’
    Authors@R:  ‘Elizabeth Eli Holmes <eli.holmes@noaa.gov>’
    
For auto-check of the email that submission comes from, I use the form that my work email uses.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Elizabeth Holmes - NOAA Federal <eli.holmes@noaa.gov>’
  
*   Suggests or Enhances not in mainstream repositories:
     marssTMB
   Availability using Additional_repositories specification:
     marssTMB   yes   https://atsa-es.r-universe.dev
     
    This is correct.
     
─  checking CRAN incoming feasibility ... [20s] NOTE
   Maintainer: 'Elizabeth Eli Holmes <eli.holmes@noaa.gov>'
   
*  New maintainer:
     Elizabeth Eli Holmes <eli.holmes@noaa.gov>
   Old maintainer(s):
     Elizabeth Holmes - NOAA Federal <eli.holmes@noaa.gov>
     
  Correct.
   
   Possibly misspelled words in DESCRIPTION:
     TMB (23:89)
     marssTMB (24:16)
   
     
## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
