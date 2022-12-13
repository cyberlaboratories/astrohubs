# this script installs mesa in the the mesahub cyberhubs application
# https://github.com/cyberlaboratories/astrohubs
#
# Two parameters can be set:
# 1: installation location (defaults to /user/mesa)
# 2: mesa version (defaults 9331)
#
echo "************** NOTE START ***************"
echo "This script will remove the link"
echo "/user/mesa/mesa (or whatever your"
echo "mesa_user_dir variable is) and set it"
echo "to the version specified below. Make"
echo "sure you avoid having this link overwritten,"
echo "for example by creating your own link manually"
echo "for example in /user/scratch14_outreach/your_user_name"
echo "and updating your MESA_DIR in yoiur .bash_aliases"
echo "************** NOTE END ***************"
echo ""
# This script has been tested on the mesahub application v1.2.6
# 
# 2018-2022, Falk Herwig, UVic
#
# define some things
mesa_user_dir=/user/mesa
# mesa_version=5329 # --> Pavel (90Msun and supermassive stars)
mesa_version=7503 # --> requested by Emma Beasor 
# mesa_version=7624 # --> Pavel (HB stars, RAWDs), Rob Farmer
# mesa_version=7184 # --> Austin (massive stars)
# mesa_version=8118 # --> Ondrea (Pop III)  
# mesa_version=8845 # --> Jacqueline (rotating models, massive stars)
# mesa_version=9331 # --> tested mesa_h5 output with
# mesa_version=10398 
# mesa_version=12115 # --> Brad RCB stars
# mesa_version=15140 # --> RCB Brad paper II


# probably nothing needs to be changed below here
mesa_source_dir=$mesa_user_dir/mesa_$mesa_version
[ ! -d $mesa_user_dir ] && mkdir $mesa_user_dir
if [ ! -d $mesa_source_dir ]  
then
    echo Installing into $mesa_user_dir
    if [ $mesa_version -lt 12115 ]
    then
        svn co -r $mesa_version svn://svn.code.sf.net/p/mesa/code/trunk $mesa_source_dir
    elif [[ $mesa_version == 15140 ]]
    then
        cd $mesa_user_dir
        wget https://zenodo.org/record/4311514/files/mesa-r15140.zip
        unzip mesa-r15140.zip
        rm mesa-r15140.zip
        mv mesa-r15140 mesa_$mesa_version
        cd mesa_$mesa_version
    else
        svn co -r $mesa_version https://subversion.assembla.com/svn/mesa^mesa/trunk $mesa_source_dir
    fi
    if [ -h $MESA_DIR ]
    then
	echo $MESA_DIR points to existing version. Remove link and set to new version $mesa_version.
	rm $MESA_DIR
    else
	echo "Make link of $MESA_DIR to version $mesa_version in /user/mesa"
    fi
    ln -s $mesa_source_dir $MESA_DIR  # $MESA_DIR is defined in .bash_aliases
    cd $mesa_source_dir/utils
    if [[ $mesa_version -lt 15140 ]]
    then
	sed -i s/"USE_PGSTAR = YES"/"USE_PGSTAR = NO"/g makefile_header
	sed -i /"LOAD_PGPLOT ="/c\ "LOAD_PGPLOT =" makefile_header
    fi
    if [[ $mesa_version =~ ^7 ]] || [[ $mesa_version =~ ^5 ]] || [[ $mesa_version =~ ^8 ]] 
    then
        sed -i '/FCbasic = -fno-range-check/s/$/ -Wno-uninitialized/' makefile_header
    fi
    cd ..
    if [[ $mesa_version =~ ^5 ]]  
    then
	cd star/private
	sed -i s/"stop'fixup'"/"stop 'fixup'"/g mod_diffusion_procs.f
	cd ../public
	sed -i s/"call do_pgstar_plots"/"\!call do_pgstar_plots"/g star_lib.f
	cd ../..
    fi
    if [[ $mesa_version == 7184 ]]
    then
	cd star/private
	sed -i s/"stop'fixup'"/"stop 'fixup'"/g mod_diffusion_procs.f
	cd ../..
    fi
    if [[ $mesa_version =~ ^9 ]] || [[ $mesa_version =~ ^8 ]] || [[ $mesa_version == 7624 ]] || [[ $mesa_version == 7503 ]]
    then
	cd star/private
	sed -i s/"stop'fixup'"/"stop 'fixup'"/g diffusion_procs.f90
	cd ../..
    fi
    ./clean
    ./install
    echo "MESA installation version $mesa_version complete in $mesa_source_dir"
else
    if [ -h $MESA_DIR ]
    then
        rm $MESA_DIR
    fi
    ln -s $mesa_source_dir $MESA_DIR  # $MESA_DIR is defined in .bash_aliases

    echo "MESA version $mesa_version already installed in $mesa_source_dir."
    echo "Force $MESA_DIR to link to" $mesa_version
fi
