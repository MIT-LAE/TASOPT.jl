# @echo off
export _calldir=$(pwd)

EMI_TOP="./opt/npss/emi-14.0.1"
export EMI_TOP

#Version String for potential version specific generated folders
NPSS_VERSION_STRING="NPSS_3.2_linux_Release_gcc5.4.0_64"
export NPSS_VERSION_STRING

. /opt/npss/npss_emi_env.sh

#export NPSS_CONFIG=nt
#export NPSS_TOP=C:\NPSS\NPSS.nt.ver32_VC14_64_LE
#export NPSS_DEV_TOP=C:\NPSS\NPSS.nt.ver32_VC14_64_LE
#export DCLOD_PATH=%NPSS_TOP%\DLMComponents\nt
#export PATH=%PATH%;%NPSS_TOP%\bin;%NPSS_TOP%\scripts
#export PATH=%EMI_TOP%\scripts\nt;%EMI_TOP%\scripts\AutoDoc;%PATH%

#Add emi scripts to path

PATH="$EMI_TOP/scripts/$NPSS_CONFIG/:$PATH"
PATH="$EMI_TOP/scripts/AutoDoc/:$PATH"
export PATH

cd $_calldir
export _calldir=''

# Run file
runnpss TEsys.run -I ./src -I ./model -I ./data