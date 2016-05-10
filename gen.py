#!/usr/bin/python2.7

out=open("./ion", "w")
with open("./ionosphere.dat", "r") as inp:
	for line in inp:
		l=list(line)
		if l[len(l)-2]=='b':
			l[len(l)-2]='1'
		else:
			l[len(l)-2]='-1'
		out.write("".join(l))
out.close()
