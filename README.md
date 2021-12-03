# domains
Script to make domains in a course grained elastic network.

# usage
Domains are defined in a domains.dat file using a simple format:
```
1
1 297 
1
298 362 
1
364 408 
1
411 428 
1
434 450 
1
501 575
```
The first number signifies the number of separate segments are in a single domain.
Then that number of lines are used to specify the _residue number_ ranges that belong to the same domain.
When run, the script will cut elastic bonds that are _between_ domains, and keep bonds _within_ domains intact.
