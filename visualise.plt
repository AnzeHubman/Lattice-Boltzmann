################################################################################
# A simple Gnuplot script to visualise output                                  #
# Two options are possible:                                                    #
#  a) for plotting velocity field use: splot 'vfield' u 1:2:3 with image       #
#  b) for plotting pressure dist. use: splot 'vfield' u 1:2:4 with image       #
#                                                                              #
# palette was found here: https://github.com/Gnuplotting/gnuplot-palettes      #
################################################################################

# line styles
set style line 1 lt 1 lc rgb '#B2182B' # red
set style line 2 lt 1 lc rgb '#D6604D' # red-orange
set style line 3 lt 1 lc rgb '#F4A582' # 
set style line 4 lt 1 lc rgb '#FDDBC7' # pale orange
set style line 5 lt 1 lc rgb '#D1E5F0' # pale blue
set style line 6 lt 1 lc rgb '#92C5DE' # 
set style line 7 lt 1 lc rgb '#4393C3' # medium blue
set style line 8 lt 1 lc rgb '#2166AC' # dark blue

# palette
set palette defined ( 0 '#B2182B',\
    	    	      1 '#D6604D',\
		      2 '#F4A582',\
		      3 '#FDDBC7',\
		      4 '#D1E5F0',\
		      5 '#92C5DE',\
		      6 '#4393C3',\
		      7 '#2166AC' )

set view map
set size ratio -1
set xrange [0:1000]
set yrange [0:201]
set xlabel 'x'
set ylabel 'y'
splot 'vfield' u 1:2:3 with image notitle
pause -1