# this script installs mesa in the the mesahub cyberhubs application
# https://github.com/cyberlaboratories/astrohubs
#
# One optional argument: installation location (defaults to /user/mesa)
#
# 2018, Falk Herwig, UVic
#
# define some things
if [ ! $1 == "" ]
then
    mesa_user_dir=$1
else
    mesa_user_dir=/user/mesa
fi
mesa_version=9331
mesa_source_dir=$mesa_user_dir/mesa_$mesa_version

# probably nothing needs to be changed below here
[ ! -d $mesa_user_dir ] && mkdir $mesa_user_dir
if [ ! -d $mesa_source_dir ] 
then
    echo Installing into $mesa_user_dir
    svn co -r $mesa_version svn://svn.code.sf.net/p/mesa/code/trunk $mesa_source_dir
    ln -s $mesa_source_dir $MESA_DIR  # $MESA_DIR is defined in .bash_aliases
    cd $mesa_source_dir/utils
    sed -i s/"USE_PGSTAR = YES"/"USE_PGSTAR = NO"/g makefile_header
    sed -i s/'LOAD_PGPLOT = `mesasdk_pgplot_link` -lz'/'LOAD_PGPLOT = '/g makefile_header
    cd ..
    ./clean
    ./install
    echo "MESA installation version $mesa_version complete in $mesa_source_dir"
else
    echo "MESA version $mesa_version already installed in $mesa_source_dir."
fi
