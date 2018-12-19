fd = 'run00.dat'
fevo = 'evolvedump.dat'
frem = 'remapdump.dat'

stats fd nooutput

do for [a = 1: int(STATS_blocks - 1)] {
  plot fd u 2:3 index a w li, fevo u 2:5 index a, frem u 2:5 index a
  pause 0.1
}
