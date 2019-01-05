alias ip='ipython --quick'
export EDITOR='emacs -nw'
alias ed="emacs -nw"

cedar() {
	ssh -Y $1@cedar.computecanada.ca
}
niagara() {
	ssh -Y $1@niagara.computecanada.ca
}
hermes() {
	ssh -Y $1@hermes.westgrid.ca
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
helix() {
	ssh $1@helix.phys.uvic.ca
}
