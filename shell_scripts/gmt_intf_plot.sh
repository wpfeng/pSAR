#!/usr/bin/env bash
#
#
#
if [ $# -lt 2 ]; then
  echo ""
  echo "$0 <intf_grd> <outps> [-gmt_proj -JM2.1i] [-gmt_ext None in default] "
  echo "                      [-wrap] [-azi -12] [-vmin -1.5] [-vmax 1.5]"
  echo "                      [-cpt jet]"
  echo "++++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "To plot an interfergoram with a bash-script using GMT tools"
  echo " -gmt_ext   -R110/112/32/33, None in default. The entire interferogram will be plotted"
  echo " -wrap      a keyword, if given, the re-wrapping will be triggered"
  echo " -vmin      -1.5 cm in default"
  echo " -vmax       1.5 cm in default"
  echo " -azi       -12 degree. if <intf_grd>.rsc exists, azi will be replaced by HEADING_DEG"
  echo " -wavelength 5.56, in default"
  echo " -gmt_proj  -JM2.1 in default"
  echo ""
  echo "Examples:"
  echo "> gmt_intf_plot.sh  geo_20240114_20240219_T56.phs.grd geo_20240114_20240219_T56.phs.ps  -vmin -20 -vmax 20 -azi -12 -gmt_proj -JM2.5i"
  echo "> "
  exit 0
fi
#
inf1=$1
outps=$2
#
gmt_proj="-JM2.1i"
#
wavelength="5.554"
gmt_ext=""
vmin="-1.5"
vmax="1.5"
wrap=0
azi="0"
cpt='jet'
PI=$(echo 'scale=50; 4*a(1)' | bc -l)
#
argv=("$@")
while [[ $k -le $# ]]; do
  case "${argv[$k]}" in
      "-gmt_proj")
	let k++
	gmt_proj=${argv[$k]}
	;;
      "-gmt_ext")
         let k++
         gmt_ext=${argv[$k]}
         ;;
      "-vmin")
          let k++
          vmin=${argv[$k]}
          ;;
      "-vmax")
          let k++
          vmax=${argv[$k]}
          ;;
      "-azi")
          let k++
          azi=${argv[$k]}
          ;;
      "-wrap")
          wrap=1
          ;;
      "-wavelength")
	  let k++
	  wavelength=${argv[$k]}
	  ;;
  esac
  let k++
done
#
rsc=`echo $inf1 | sed 's/\.grd$/\.rsc /g'`
#
if [ -e ${rsc} ]; then
  heading=`grep HEADING_DEG $rsc|awk '{print $2}'`
  if [ $azi == "0" ]; then
	 azi=$heading 

  fi
  #
  if [ $wavelength == "0" ]; then
	wavelength=`grep WAVELENGTH $rsc|awk '{print $2*100}'`
	#
  fi
fi
#
#
gmt begin $outps eps,pdf,png
  #
  gmt gmtset MAP_ANNOT_OBLIQUE 32
  gmt gmtset FONT_LABEL 8p,Helvetica,black
  #
  inv=`echo $vmin $vmax|awk '{print ($2-$1)/100}'`
  #
  echo " gmt makecpt -C$cpt -T$vmin/$vmax/$inv -D"
  gmt makecpt -C$cpt -T$vmin/$vmax/$inv -D
  #
  basename=$(basename $inf1 .grd)
  #
  if [ $wrap == 1 ]; then
     #
     ingrd="$basename.rewrap.grd"
     #
     echo " pGRD_rewrap.py $inf1 $ingrd -wavelength $wavelength -vmin $vmin -vmax $vmax"
     if [ ! -e $ingrd ]; then
        pGRD_rewrap.py $inf1 $ingrd -wavelength $wavelength -vmin $vmin -vmax $vmax
     fi 
     #
     #
  else
     ingrd="$basename.los.grd"
     #
     echo " gmt grdmath $inf1 1 NaN IFELSE  = mask.grd"
     gmt grdmath $inf1 1 NaN IFELSE  = mask.grd
     #
     echo " gmt grdmath $inf1 -1 MUL $wavelength MUL 4 DIV $PI DIV mask.grd MUL = $ingrd"
     gmt grdmath $inf1 -1 MUL $wavelength MUL 4 DIV $PI DIV mask.grd MUL  = $ingrd
  fi
  #
  echo " gmt grdimage $ingrd $gmt_proj $gmt_ext -Bxaf -Byaf -BWSne"
  gmt grdimage $ingrd $gmt_proj $gmt_ext -Bxaf -Byaf -BWSne
  #
  gmt set MAP_TICK_PEN_PRIMARY 0.25p,black
  gmt set FONT_ANNOT_PRIMARY 8p FONT_LABEL 10p
  gmt set MAP_FRAME_PEN thinner,black
  echo " gmt colorbar -DjBR+w0.75i/0.125c+o0.3c/0.3i+h+ml -BWSEN -Bxaf+l'LOS(cm)'"
  gmt colorbar -DjBR+w0.75i/0.125c+o0.3c/0.3i+h+ml -BWSEN -Bxaf+l"LOS(cm)" #-G-3.14/3.14
  #
  gmt set MAP_FRAME_PEN thicker,black
  #
  gmt inset begin -R0/1.1/0/1.1 -JX2c/2c -DjBL+o0.01c/0.01c
    # start 0.75,0
    abs_azi=`echo $azi |awk '{if ($1 < 0) print -$azi; else $azi}'` 
    flag=`echo $abs_azi |awk '{if ($1 > 150) print 1; else print 0}'`
    #
    if [ $flag -eq 1 ]; then
       x0="0.9"
       y0="0.9"
    else
       x0="0.75"
       y0="0"
    fi
    #
    #
    vx=`echo $x0 $azi |awk '{print sin($2/180*3.1415926)}'`
    vy=`echo $y0 $azi |awk '{print cos($2/180*3.1415926)}'`
    echo $x0 $y0 $vx $vy 0 0 0.001 > arrow_flight.gmt
    gmt velo arrow_flight.gmt -Se1.5c/0.0+f0 -A0.5c+e+p0.0p,black -W3p,black -Gblack
    #
    #
    x1=`echo $x0 $azi |awk '{print sin($2/180*3.1415926)*0.35+$1}'`
    y1=`echo $y0 $azi |awk '{print cos($2/180*3.1415926)*0.35+$1}'`
    #
    vx=`echo $x0 $azi |awk '{print 0.5*sin(($2+90)/180*3.1415926)}'`
    vy=`echo $x0 $azi |awk '{print 0.5*cos(($2+90)/180*3.1415926)}'`
    echo $x1 $y1 $vx $vy 0 0 0.001 > arrow_range.gmt
    gmt velo arrow_range.gmt -Se1.5c/0.0+f0 -A0.5c+e+p0.0p,black -W2p,black -Gblack
    #
  gmt inset end
  #
gmt end

