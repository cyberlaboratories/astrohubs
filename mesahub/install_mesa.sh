# this script installs mesa in the the mesahub cyberhubs application
# https://github.com/cyberlaboratories/astrohubs
#
# Two parameters can be set:
# 1: installation location (defaults to /user/mesa)
# 2: mesa version (defaults 9331)
#
# 2018, Falk Herwig, UVic
#
# define some things
mesa_user_dir=/user/mesa
#mesa_version=7624 # --> Pavel (HB stars, RAWDs), Rob Farmer
#mesa_version=7184 # --> Austin (massive stars)
mesa_version=8118 # --> Ondrea (Pop III)
#mesa_version=8845 # --> Jacqueline (rotating models, massive stars in the future)
#mesa_version=9331 # --> tested mesa_h5 output with
#mesa_version=10398 
# note: there may be some cyberhub images still around that have problems with -lz libaray, please report to fherwig@uvic.ca

# probably nothing needs to be changed below here
mesa_source_dir=$mesa_user_dir/mesa_$mesa_version
[ ! -d $mesa_user_dir ] && mkdir $mesa_user_dir
if [ ! -d $mesa_source_dir ] 
then
    echo Installing into $mesa_user_dir
    svn co -r $mesa_version svn://svn.code.sf.net/p/mesa/code/trunk $mesa_source_dir
    ln -s $mesa_source_dir $MESA_DIR  # $MESA_DIR is defined in .bash_aliases
    cd $mesa_source_dir/utils
    sed -i s/"USE_PGSTAR = YES"/"USE_PGSTAR = NO"/g makefile_header
    sed -i /"LOAD_PGPLOT ="/c\ "LOAD_PGPLOT =" makefile_header
    if [[ $mesa_version =~ ^7 ]]
    then
        sed -i '/FCbasic = -fno-range-check/s/$/ -Wno-uninitialized/' makefile_header
    fi
    cd ..
    ./clean
    ./install
    echo "MESA installation version $mesa_version complete in $mesa_source_dir"
else
    if [ -h $MESA_DIR ]
    then
	rm $MESA_DIR
	ln -s $mesa_source_dir $MESA_DIR  # $MESA_DIR is defined in .bash_aliases
    fi
    
    echo "MESA version $mesa_version already installed in $mesa_source_dir."
    echo "$MESA_DIR points to" $mesa_version
fi
