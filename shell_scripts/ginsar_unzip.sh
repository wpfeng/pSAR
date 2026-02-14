#!/usr/bin/env bash
#
#
# Developed by Feng, W.P., @NRCan, 2015-11-06
#  rewrite it for unzip zips 
#
#############################################
#

if [ $# -lt 1 ]; then
   cat << EOF
   
   $0 <obj_dir> [searchstring, *.zip in default]
   +++++++++++++++++++++++++++++++++++++++++++++++++++++
   Developed by Wanpeng Feng, @NRCan, 2015-11-06
   Updated by Wanpeng Feng, @NRCan, 2016-11-28

EOF
   exit -1
fi
searchstr="*.zip"
if [ $# -lt 1 ]; then
   objdir=$pwd
else
   objdir=$1
fi
if [ $# -ge 2 ]; then
   searchstr=$2
fi
#
#
echo " ******************************************************"
echo " + gInSAR: Uncompress ZIP(s) in the current folder"
echo " + gInSAR: $objdir and $searchstr"
echo " ******************************************************"
echo ""
#
for czip in `ls $searchstr -f`; do
  #
  cfile=`basename $czip .zip`
  #
  echo "  $objdir/$cfile"
  if [ ! -d $objdir/$cfile ]; then
     #
     mkdir $objdir/$cfile
     unzip -n $czip -d $objdir/$cfile
     #
  fi 
  #
done

exit 0
