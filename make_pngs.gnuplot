#!/usr/bin/env gnuplot
# USAGE EXAMPLE:
#   gnuplot -e "start=0; end=300; step=10" make_pngs.gnuplot

# default values if not provided
if (!exists("start"))   start = 0
if (!exists("end"))     end = 9900
if (!exists("step"))    step = 10

unset key
set cbrange [0:0.00015]
set terminal pngcairo size 800,600
# Generate individual PNG frames
do for [i=start:end:step] {
    filename = sprintf("out_%d.data", i)
    if (system(sprintf("test -f %s", filename))) {
        print sprintf("Skipping missing file: %s", filename)
        continue
    }
    outname = sprintf("frames/out_%d.png", i)
    system("mkdir -p frames")
    set output outname
    print sprintf("Saving %s", outname)
    plot filename with image
    unset output
}