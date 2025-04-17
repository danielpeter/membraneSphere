#!/usr/bin/env python

from Numeric import *
import Gnuplot, Gnuplot.funcutils

# start gnuplot	
#g = Gnuplot.Gnuplot(persist=1)
g = Gnuplot.Gnuplot(debug=0)

#view options
g('set terminal x11')
#g('set terminal jpeg')
g('unset colorbox')
g('unset key')
g('unset border')
g('set xtics nomirror')
g('set ytics nomirror')
g('unset xtics')
g('unset ytics')
g('unset ztics')

g('set cbrange [-100:100]')
g('set xrange[0:8000]')  #[0:2000]
g('set yrange [0:8000]') #[0:2000]
g('set zrange [-4000:4000]') 
g('set view 147,123,1.3,1.0')

#g.title('simulation')
	
for n in arange(0,300,10):   #(0,-90,-1) #delta.dr8
 #file
 time=str(n)
 time='%(#)04d' % {"#":n}  
 filename='simulation.'+time+'.dat' 

 #show
 #g('set output "Data/simulation/movie/'+filename+'.jpg" ')

 g.splot('"'+filename+'" u ($1*(6000+$4)):($2*(6000+$4)):($3*(6000+$4)):4 w p palette')
 print time
 
 # print output 
 #g('set terminal pdf')
 #g('set output "Data/simulation/movie/'+filename+'.pdf" ')
 #g('replot')
 #g('set terminal x11')
 #g('set output')
 
 #g('pause 0.1')
