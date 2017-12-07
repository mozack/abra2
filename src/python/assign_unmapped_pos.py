import sys

for line in sys.stdin:
	line = line.rstrip()
	if line.startswith('@'):
		print line
	else:
		fields = line.split()
		flag = int(fields[1])
		# 4 = read unmapped.  8 = mate unmapped
		if flag & 0xC == 4 and fields[2] == '*':
			# This read is unmapped and mate is mapped
			chrom = fields[6]
			pos = fields[7]
			fields[2] = chrom
			fields[3] = pos
			print '\t'.join(fields)
		elif flag & 0xC == 8 and fields[6] == '*':
			# This read is mapped and mate is unmapped
			chrom = fields[2]
			pos = fields[3]
			fields[6] = chrom
			fields[7] = pos
			print '\t'.join(fields)
		else:
			print line