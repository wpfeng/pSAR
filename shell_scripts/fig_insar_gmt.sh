#!/usr/bin/env bash
#
#
if [ $# -lt 1 ]; then
   echo 
   echo "  $0 <input.grd>  <outps> [-azi <azi>] "
   echo "                  [-gext minlon/maxlon/minlat/maxlat] "
   echo "                  [-J M] [-size 10] [-cpt jet] [-cpt_rng -1/1/0.01] "
   echo "                  [-epi lon,lat] [-wave 5.54 -unit cm] [-scale 20] [-wrap 3.1415926]"
   echo "  +++++++++++++++++++++++++++++++++++++++++++++++"
   echo "  Plot the (GRD)files in GMT."
   echo "  "
   echo "  Example:"
   echo "   > fig_insarfitting_gmt.sh geo_20240114_20240126_T56.phs.grd geo_20240114_20240126_T56.ps -wrap 5 -cpt_rng -5/5/0.01 -azi -12 -scale 8"
   echo "   " 
   exit 0
fi
#
# an internal function to plot an interferogram
#
img_subplot(){
	#
  #
  local ingrd=$1
  local cpt_phs=$2
  local cpt_phs_rng=$3
  local gext=$4
  local peojection=$5
  local size=$6
  local label=$7
  local azi=$8
  local xshift=${9}
  local yshift=${10}
  local Bnswe=${11}
  local epi=${12}
  #
  color_length=`echo $size|awk '{print $1*4/10}'`
  color_width=`echo $size|awk '{print $1*0.5/10}'`
  inset_width=`echo $size|awk '{print $1*2/10}'`
  font_size=`echo $size|awk '{print $1*18/10}'`
  font_size1=`echo $size|awk '{print $1*18/10-2}'`
  font_size2=`echo $size|awk '{print $1*18/10-4}'`

  epi_size=`echo $size|awk '{print $1*1/40}'`
  #

  gmt set MAP_ANNOT_OBLIQUE 34
  # gmt set MAP_ANNOT_ORTHO z
  gmt set MAP_FRAME_TYPE plain
  gmt set FONT Helvetica,black
  gmt set FONT_ANNOT_PRIMARY 8p FONT_LABEL 10p
  gmt set MAP_TICK_PEN_PRIMARY 0.125p,black
  gmt set MAP_ANNOT_OFFSET_PRIMARY 0.25p
  gmt set FONT_ANNOT_PRIMARY Times-Roman
  gmt set FONT_LABEL Times-Roman
  gmt set FONT_TITLE Times-Roman
  #
  gmt makecpt -C$cpt_phs -T$cpt_phs_rng -D -Z
  #
  #
  # echo " gmt grdimage -J"$projection$size"c $ingrd -Q -X$xshift -Y$yshift  -Bxaf -Byaf $Bnswe"
  gmt grdimage $ingrd -J"$projection$size"c -Q -X$xshift -Y$yshift -Bxaf -Byaf $Bnswe
  #gmt plot $fault_model -W2.0p,black

  echo $epi| gmt plot -Sa${epi_size}i -W0.8p,black,solid -Gred
  #
  #echo "$label" | gmt text -F+a0+f${font_size}p,Helvetica,white+cTL -Gblack -Dj0.05c/0.1c
  #
  gmt set MAP_TICK_PEN_PRIMARY 0.5p,black
  gmt set FONT_ANNOT_PRIMARY 8p FONT_LABEL 10p
  #
  gmt colorbar -DjBL+w${color_length}c/${color_width}c+o0.2c/0.60c+h+ml -Bxaf+l"LOS ($unit)" --FONT_ANNOT_PRIMARY=${font_size2}p --FONT_LABEL=${font_size1}p --MAP_FRAME_PEN=1p
  #
  gmt basemap -LjBL+w$scale+u+o0.65c/1.75c
  #
  gmt set MAP_TICK_PEN_PRIMARY 0.8p,black
  gmt set FONT_ANNOT_PRIMARY 14p FONT_LABEL 15p
  #
  ##arrow
  gmt inset begin -R0/2/0/2 -JX${inset_width}c/${inset_width}c -DjBR+o0.1c/0.1c
    #
    abs_azi=`echo $azi |awk '{if ($1 < 0) print -$azi; else $azi}'`
    flag=`echo $abs_azi |awk '{if ($1 > 150) print 1; else print 0}'`
    if [ $flag -eq 1 ]; then
       x0="1.75"
       y0="1.95"
    else
       x0="1.25"
       y0="0.15"
    fi
    #
    vx=`echo $x0 $azi |awk '{print sin($2/180*3.1415926)}'`
    vy=`echo $y0 $azi |awk '{print cos($2/180*3.1415926)}'`
    echo $x0 $y0 $vx $vy 0 0 0.001 > arrow_flight.gmt
    gmt velo arrow_flight.gmt -Se1.0c/0.0/0.8+f0 -A0.35c+e+p0.0p,black -W2p,black -Gblack
    #
    x1=`echo $x0 $azi |awk '{print sin($2/180*3.1415926)*0.65+$1}'`
    y1=`echo $y0 $azi |awk '{print cos($2/180*3.1415926)*0.65+$1}'`
    #
    vx=`echo $x0 $azi |awk '{print 0.5*sin(($2+90)/180*3.1415926)}'`
    vy=`echo $x0 $azi |awk '{print 0.5*cos(($2+90)/180*3.1415926)}'`
    echo $x1 $y1 $vx $vy 0 0 0.001 > arrow_range.gmt
    gmt velo arrow_range.gmt -Se1.0c/0.0/0.8+f0 -A0.35c+e+p0.0p,brown -W1.5p,brown -Gbrown
    #
    #
  gmt inset end
}
#
raw_grd=$1
outps=$2
cpt_phs="jet"
cpt_phs_rng="-1.5/1.5/0.01"
#
epi="87.45 28.5"
gext=$raw_grd
projection="M"
size=6
wavelength=5.546
minv=-3
maxv=3
unit="cm"
azi="-169"
scale=40
wrap=None
#
argv=("$@")
k=0
while [[ $k -le $# ]]; do
  case "${argv[$k]}" in
      "-azi")
	      let k++
          azi=${argv[$k]}
          ;;
      "-gext")
          let k++
          gext=${argv[$k]}
          ;;
      "-cpt")
          let k++
          cpt_phs=${argv[$k]}
          ;;
      "-cpt_rng")
          let k++
          cpt_phs_rng=${argv[$k]}
          ;;
      "-fault")
          let k++
          fault_model=${argv[$k]}
          ;;
      "-epi")
          let k++
          epi="${argv[$k]} ${argv[$k+1]}"
          let k+=1
          ;;
      "-wave")
          let k++
          wavelength=${argv[k]}
          ;;
      "-unit")
          let k++
          unit=${argv[k]}
          ;;
      "-minv")
          let k++
          minv=${argv[k]}
          ;;
      "-maxv")
          let k++
          maxv=${argv[k]}
          ;;
      "-J")
          let k++
          projection=${argv[k]}
          ;;
      "-size")
          let k++
          size=${argv[k]}
          ;;
      "-scale")
          let k++
          scale=${argv[k]}
          ;;
      "-wrap")
          let k++
          wrap=${argv[k]}
          ;;
  esac
  let k++
done
#
# convert from radian to light-of-sight displacements...
#
gmt grdmath -$wavelength 4 DIV PI DIV $raw_grd MUL = tmp_raw_los.grd
#
if [ "$wrap" != "None" ]; then
    # 
    # gmt grdmath input.nc a ADD 2 a MUL FMOD a SUB = output.nc
    cpt_phs_rng="-$wrap/$wrap/0.001"
    #
    gmt grdmath tmp_raw_los.grd $wrap ADD 2 $wrap MUL MOD $wrap SUB = tmp_raw_los.grd
else
    echo "We don't need rewrap the intfs."
fi
###########
gmt begin $outps png,pdf
  #
  ############################# Plot RAW_file.grd
  #
  echo " img_subplot tmp_raw_los.grd ${cpt_phs} $cpt_phs_rng $gext $projection $size (a) $azi 1i 1i -BnSWe [$epi]"
  img_subplot tmp_raw_los.grd "$cpt_phs" "$cpt_phs_rng" "$gext" "$projection" "$size" "(a)" "$azi" "1i" "1i" -BnSWe "$epi"
  #
  rm tmp_*.grd
gmt end
