# Command line options
```
Usage: fxv [-options] [output] input
   input              if only input then extract
  -v                  show verbose information
  -0                  store
  -1                  compress file
  -h                  extended help
  -2                  tune on file, output file is not created
  -t<n>               n is number of threads, default=1
  -p<component><idx>  select tunable component:
                      a - SMC    b - APM1    c - DS     d - AVG
                      e - SCM    f - RCM     g - CM     h - MX
                      i - ST     j - DHS     k - SM     l - SK
                      m - APM2   n - ERR     o - UAS    p - LMX
                      q - STA    r - BYT
                      idx - select component index
  -o<n>               n specifies percentage of tune, default=100
  -r<n>               number of tune runs, default=25
  -m<n>               minimum tune improvment in bytes, default=2
  -z<n>               select stream for tuning, use after -2 option, default=all
  -f<n>               tune all parameters, use n paramters max, default=all
  -bc                 bc - enable bounds check at compile, dafault=false
  -br                 br - enable bounds check at runtime, dafault=false
  -k                  k - disable radius in tune, dafault=true
  -j                  j - do x86 JIT, dafault=false
  -i                  i - show cfg component info, default=false
  -w                  i - do not store stream model in archive
  -s<n>               s - seed for tune, default=random
  -c<file>            c - use config file. dafault=config.pxv
  -d dir1/input       extract to dir1
  -d dir1/input dir2  extract to dir2
  -l input            list archive
```
### option -p
For example, selecting the BYT indices of the component: -pr0011

Tune only BYT component 2 (third value is 1) and 3 (fourth value is 1). BYT 0 (first value is 0) and 1 (second value is 0) are disabled.

### logging
Program creates a pxv.log file in the same directory. After compression finishes, the final result is written at the end of the file.
Example line from log:
```
enwik9 1000000000 -> 641256112 (5.1300 bpc) in 555.38 sec (1758.381 KB/sec), 7432 Kb  Tue Nov 26 19:34:15
```
