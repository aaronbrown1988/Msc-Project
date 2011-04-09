set terminal gif
set output "OUTNAME.gif"
set view 0, 0, 1,1
set title "NAME"
set pm3d map
splot "map-NAME.tsv" 
