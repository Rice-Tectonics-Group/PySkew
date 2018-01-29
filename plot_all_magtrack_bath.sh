#!/bin/bash

#In file defs
infile=$1
grdfile=ETOPO1_Bed_g_gmt4.grd
outdir="`dirname $infile`/bathymetry"

if [ ! -d "$outdir" ]; then
    mkdir $outdir
fi

#read deskew file
while read i; do
    if [ "${i:0:9}" != "comp_name" ]; then
        echo $i | awk '{if (substr($1,1,1) != "#") printf("%s\t%s\t%s\t%s/%s.deskewed\t%s\n", $11, $6, $7, $9, $1, $8-180)}' >> tmp.txt
    fi
done < $infile

#Flags and parameters
proj=-Jm0.5c
cpt=-Cnew.cpt
cptfile=new.cpt
zscale=-Z300
sitedis=10.0
sitesize=.15 #in degrees unless c,i,p appended at end to mean cm, inch, point respectivly
sitecolor=limegreen #rgb, hex, or name
sitealpha=15 #higher number=more transparent

#Script start

 makecpt -Cpolar -T-6000/0/400 > $cptfile

while read t slat slon f azi; do
    if [ $t == "aero" ]; then
        color="black"
    else
        color="darkslategrey"
    fi

    #outfile path in data file option
#    plot="`dirname $f` /bath.ps"

    #outfile path in single bath file
    plot="$outdir/`basename $f`.ps"
    plotpng="`dirname $plot`/`basename $plot .ps`.png"
    plotpdf="`dirname $plot`/`basename $plot .ps`.pdf"

    #Map boundary defs
    y1=$(echo "$slat-$sitedis" | bc)
    y2=$(echo "$slat+$sitedis" | bc)
    x1=$(echo "$slon-$sitedis" | bc)
    x2=$(echo "$slon+$sitedis" | bc)
    region=-R$x1/$x2/$y1/$y2

     grdimage $grdfile $cpt $proj $region -K > $plot

     pscoast $region $proj -W1p,black -Ggrey -B5 -K -O >> ${plot}

     psscale -Dn1.2/0+w10/1 -O -K $region $proj $cpt >> ${plot}

    echo "$slon $slat" |  psxy $region $proj -O -K -Sc$sitesize -G$sitecolor -t$sitealpha >> ${plot}

     pswiggle $f $region $proj $zscale -O -I$azi -G-$color -Wthin,$color -Tthin,$color,- -t50 >> ${plot}

    ps2pdf $plot $plotpdf
#    convert $plot $plotpng #possible png option???

    echo "finishing plotting $plotpdf"

done < tmp.txt

rm tmp.txt
rm new.cpt
rm $outdir/*.ps


