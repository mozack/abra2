import re
import sys

# Filter inserts at alignment start or end that are adjacent to intron
# Filter deletions adjacent to intron
#p = re.compile('^[0-9]*I[0-9]*N.*|^[0-9]*S[0-9]*I[0-9]*N.*|.*N[0-9]*I$|.*N[0-9]*I[0-9]*S$|.*N[0-9]*[D].*|.*D[0-9]*N.*')

# More aggressive filter.  Do not require inserts to be at read start or end
# Filter inserts adjacent to intron
# Filter deletions adjacent to intron
p = re.compile('.*N[0-9]*[ID].*|.*[ID][0-9]*N.*')


for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('@'):
        print line
        continue
    
    shouldFilter = False
         
    fields = line.split()
    cigar = fields[5]
    
    m = p.match(cigar)
    if m:
        shouldFilter = True
    
    idx = 11
    while idx < len(fields):
        tag = fields[idx]
        if tag.startswith('MC:Z'):
            m = p.match(tag)
            if m:
                shouldFilter = True
                break
        idx += 1
     

    if not shouldFilter:
        print line

