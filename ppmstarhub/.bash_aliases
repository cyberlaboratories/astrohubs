export PS1='\[\e]2;\W\a\]${PWD#"${PWD%/*/*}/"}\$ '
export EDITOR='emacs -nw'

cedar() {
        ssh -Y $1@cedar.computecanada.ca
}
niagara() {
        ssh -Y $1@niagara.computecanada.ca
}
scandium() {
        ssh -Y $1@206.12.89.164
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
alias astrohubc='ssh  centos@206.12.91.2'
alias astrohub='ssh  centos@206-12-89-142.cloud.computecanada.ca'

alias ip='ipython --quick'
alias ed="emacs -nw"
alias ec="emacsclient"

