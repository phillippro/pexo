#!/bin/sh

#################################################################################
#
# To run PEXO via this script, you need to set an environment variable $PEXODIR
# to a path to the PEXO repository.
# It is also recommended to create an alias for this script to run it from
# anywhere in the terminal.
#
# e.g.:
#
# BASH, add to ~/.bashrc or ~/.bash_profile:
#     export PEXODIR=/example/path/to/pexo
#     alias pexo="/example/path/to/pexo/pexo.sh"
#
# TCSH, add to ~/.tcshrc:
#     setenv PEXODIR /example/path/to/pexo
#     alias pexo /example/path/to/pexo/pexo.sh
#
#################################################################################

usage () {
    # usage information
    echo "$(basename "$0") [ARGUMENTS] -- a bash wrapper for PEXO software. See documentation for full reference (https://github.com/phillippro/pexo)."
    echo
    echo "Arguments:"
    echo
    echo "-m, --mode         PEXO mode: emulate or fit [optional, default='fit']"
    echo
    echo "-c, --component    PEXO model component: timing (T astrometry (A radial velocity (R) and their combinations [optional, default='TAR']"
    echo
    echo "-i, --ins          Instrument or observatory [mandatory for emulation mode; default=NA]"
    echo
    echo "-P, --par          Parameter file: parameters for astrometry and observatory [default: automatically obtained from simbad if not find parameter file]"
    echo
    echo "-N, --Niter        Length of MCMC [optional, default=1e3]"
    echo
    echo "-C, --Companion    Companion data directory: directory with timing, RV or astrometry data files [optional, default='']"
    echo
    echo "-g, --geometry     Geometric orbit or relativistic orbit [optional, default=TRUE]"
    echo
    echo "-n, --ncore        Number of cores [optional, default=4]"
    echo
    echo "-t, --time         Timing file: epochs or times could be in 1-part or 2-part JD[UTC] format [mandatory if mode=emulate, default='2450000 2460000 100']"
    echo
    echo "-p, --primary      Primary star name [mandatory]"
    echo
    echo "-s, --secondary    Secondary star name [optional, default=NULL]"
    echo
    echo "-M, --mass         Mass of primary in unit of solar mass [optional, default=1]"
    echo
    echo "-d, --data         Data directory: directory with timing, RV or astrometry data files [mandatory if mode=fit, default=NULL]"
    echo
    echo "-v, --var          Output variables [optional, default='BJDtcb BJDtdb RvTot ZB'],"
    echo
    echo "-o, --out          Output file name: relative or absolute path [optional, default=out_pexo.txt]"
    echo
    echo "-f, --figure       Output figure: FALSE or TRUE [optional, default=FALSE]"
    echo
    echo "-V, --verbose      Verbose: FALSE or TRUE [optional, default=FALSE]"
}

# parse arguments, consistent with pexo.R
while [ $# -gt 0 ]
do
    key="$1"
    case $key in
        -m|--mode) #####
        mode="$2"
            shift
            shift
            ;;
        -c|--component) #####
        component="$2"
            shift
            shift
            ;;
        -i|--ins) #####
        ins="$2"
            shift
            shift
            ;;
        -P|--par) #####
        ins="$2"
            shift
            shift
            ;;
        -N|--Niter) ####
        Niter="$2"
            shift
            shift
            ;;
        -C|--Companion) ####
        companion="$2"
            shift
            shift
            ;;
        -g|--geometry) ####
        geometry="$2"
            shift
            shift
            ;;
        -n|--ncore) ####
        ncore="$2"
            shift
            shift
            ;;
        -t|--time) ####
        time="$2"
            shift
            shift
            ;;
        -p|--primary) ####
        primary="$2"
            shift
            shift
            ;;
        -s|--secondary) ####
        secondary="$2"
            shift
            shift
            ;;
        -M|--mass) ####
        mass="$2"
            shift
            shift
            ;;
        -d|--data) ####
        data="$2"
            shift
            shift
            ;;
        -v|--var)
        var="$2"
            shift
            shift
            ;;
        -o|--out)
        out="$2"
            shift
            shift
            ;;
        -f|--figure)
        figure="$2"
            shift
            shift
            ;;
        -V|--verbose)
        verbose="$2"
            shift
            shift
            ;;
        -h|--help)
        usage
            exit 0
            ;;
        *)    # unknown option
            echo "Unknown option: $1"
            echo "Run 'pexo.sh --help' to see available arguments."
            exit 1 # past argument
            ;;
    esac
done

# check if R is installed
rscript=$(which Rscript)
if [ -z "$rscript" ]
then
    echo "Error: \$Rscript is not installed. Exiting."
    exit 1
else
    echo "Found Rscript in $rscript"
fi

# check if the PEXODIR variable is set
if [ -z "$PEXODIR" ]
then
    echo "Error: \$PEXODIR is not set. It should point to PEXO repository. Exiting."
    exit 1
else
    echo "Found PEXO in $PEXODIR"
fi

# save current path
original_path="$PWD"

# convert input paths to absolute
#time=$(realpath $time 2>/dev/null)
#par=$(realpath $par 2>/dev/null)
#out=$(realpath $out 2>/dev/null)

# go to pexo code
cd "$PEXODIR"
cd "code"

# base command
command="$rscript pexo.R"
# command="Rscript pexo.R"

# add parameters (could do with a function here, but bash functions and strings don't mix well)
if [ ! -z "$mode" ]
then
    command="$command -m $mode"
fi

if [ ! -z "$component" ]
then
    command="$command -c $component"
fi

if [ ! -z "$ins" ]
then
    command="$command -i $ins"
fi

if [ ! -z "$par" ]
then
    command="$command -P $par"
fi

if [ ! -z "$Niter" ]
then
    command="$command -N $Niter"
fi

if [ ! -z "$companion" ]
then
    command="$command -C $companion"
fi

if [ ! -z "$geometry" ]
then
    command="$command -g $geometry"
fi

if [ ! -z "$ncore" ]
then
    command="$command -n $ncore"
fi

if [ ! -z "$time" ]
then
    command="$command -t $time"
fi

if [ ! -z "$primary" ]
then
    command="$command -p $primary"
fi

if [ ! -z "$secondary" ]
then
    command="$command -s $secondary"
fi

if [ ! -z "$mass" ]
then
    command="$command -M $mass"
fi

if [ ! -z "$data" ]
then
    command="$command -d $data"
fi

if [ ! -z "$var" ]
then
    command="$command -v $var"
fi

if [ ! -z "$out" ]
then
    command="$command -o $out"
fi

if [ ! -z "$figure" ]
then
    command="$command -f $figure"
fi

if [ ! -z "$verbose" ]
then
    command="$command -V $verbose"
fi

# run PEXO
echo "---------------------------------"
eval $command

# go back to the original directory
cd $original_path
