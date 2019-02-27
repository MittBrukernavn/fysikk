# encoding: utf-8
import matplotlib.pyplot as plt
import numpy
from math import ceil

def eulerNewtons2ndLaw(T, d, s_0, v_0):
	N = ceil(T/d)
	t = numpy.zeros(N+1)
	s = numpy.zeros(N+1)
	v = numpy.zeros(N+1)	# her er det strengt tatt tilstrekkelig
	t[0] = 0				# å lagre siste verdi av v,
	s[0] = s_0				# med mindre man vil plotte den også
	v[0] = v_0
	for n in range(N):
		k = 0.5			# eks. med fjær
		a = -k*s[n]		# formel for akselerasjonen
		v[n+1] = v[n] + d*a
		s[n+1] = s[n] + d*v[n+1]	# usikker på om de egentlig 
		t[n+1] = t[n] + d			# vil ha v[n] der
	plt.plot(t, s)	# tid langs x-akse, strekning langs y
	plt.xlabel('t')
	plt.ylabel('s')
	plt.plot(t, v)
	plt.show()
	

eulerNewtons2ndLaw(100, 0.1, 1, 0)
