# rsync_monitor
Prototype code that monitors the images being transferred via rsync
and ingested into Gen2 Butler repositories.

Writes an html page for the camera + src site + dayobs.
Adds a row to a summary html page for the camera + src site.
Keeps monitoring information in a database to make updates faster.


## Requires
* Python3
* configargparse
* sqlite3
* Read access to rsync'd files, rsync logs, Gen2 repo, src site logs

## Examples

Using configargparse allows configuration values to be set on the 
command line, a yaml config file (see config dir for examples)
as well as environment variables (start with prefix RSYMON_)

The example config files assume that the following environment variables
are set: `RSYMON_DB_ROOT` and `RSYMON_HTML_ROOT`
For example:
```
    export RSYMON_DB_ROOT=${HOME}/public_html/rsync_monitor/db
    export RSYMON_HTML_ROOT=${HOME}/public_html/rsync_monitor
```

Example command lines:
* Creating new sqlite3 file (--init_db):
  ```
  ./monitor_rsync.py --init_db -c config/auxTel_L1Archiver_monitor.yaml 2019-09-30
  ```

* Running over a range of dates:
  ```
  ./monitor_rsync.py -c config/auxTel_L1Archiver_monitor.yaml 2019-09-01 2019-10-01
  ```

## Notes
If you see an error like the following, the config file is probably 
missing quotes around string value that contains curly braces.
```
    monitor_rsync.py: error: Couldn't parse config file: while parsing a block mapping
      in "config/auxTel_DAQ_monitor.yaml", line 1, column 1
   expected <block end>, but found '<scalar>'
      in "config/auxTel_DAQ_monitor.yaml", line 11, column 21
```
