fd = 'run00.dat'

stats fd nooutput

set xrange [0:1]
set yrange [0:1]

do for [a = 1: int(STATS_blocks - 1)] {
  plot fd u 2:3 index a w li
  pause 0.04
}
