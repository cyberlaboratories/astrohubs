alias ed="emacs -nw"
alias ppm2png='(names=`ls *ppm`; for name in $names; do convert $name ${name:0:21}png; rm $name; done)'
cedar() {
	ssh -Y $1@cedar.computecanada.ca
}
bluewaters() {
	ssh -Y  $1@h2ologin.ncsa.illinois.edu
}
orcinus() {
	ssh -Y $1@seawolf3.westgrid.ca
}
aquila() {
	ssh -Y $1@aquila.phys.uvic.ca
}
frodo() {
	ssh $1@frodo.lcse.umn.edu
}
