n = []
bins = 0
nbins = 20
increment = 4300/nbins
for i in range(5001):
	if i>=700 and (i-700)%increment==0:
		n.append(i)
		bins += 1
print n
print bins
