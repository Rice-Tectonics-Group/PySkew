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

gmt makecpt -Cpolar -T-5000/-2000/150 > $cptfile

while read t slat slon f azi; do
    if [ $t == "aero" ]; then
        color="black"
        alpha=70
        color2="maroon1"
        alpha2=50
        f2=`echo "$f" | sed -r 's/[E]+/V/g'`
        if [ $f == $f2 ]; then
            continue
        fi
    else
        color="black"
        alpha=50
    fi

    #outfile path in data file option
#    plot="`dirname $f` /bath.ps"

    #outfile path in single bath file
    plot="$outdir/`basename $f`.ps"
    title=`basename $f`
    plotpng="`dirname $plot`/`basename $plot .ps`.png"
    plotpdf="`dirname $plot`/`basename $plot .ps`.pdf"

    #Map boundary defs
    y1=$(echo "$slat-$sitedis" | bc)
    y2=$(echo "$slat+$sitedis" | bc)
    x1=$(echo "$slon-$sitedis" | bc)
    x2=$(echo "$slon+$sitedis" | bc)
    region=-R$x1/$x2/$y1/$y2

    gmt grdimage $grdfile $cpt $proj $region -K > $plot

    gmt pscoast $region $proj -W1p,black -Ggrey -B5 -K -O >> ${plot}

    gmt psscale -Dn1.2/0+w10/1 -O -K $region $proj $cpt >> ${plot}

    echo "$slon $slat" | gmt psxy $region $proj -O -K -Sc$sitesize -G$sitecolor -t$sitealpha >> ${plot}

    if [ $t == "aero" ]; then
        gmt pswiggle $f2 $region $proj $zscale -O -K -I$azi -G-$color2 -Wthin,$color -Tthin,$color2,- -t$alpha2 -B+t$title >> ${plot}
    fi

    gmt pswiggle $f $region $proj $zscale -O -I$azi -G-$color -Wthin,$color -Tthin,$color,- -t$alpha -B+t$title >> ${plot}

    ps2pdf $plot $plotpdf
#    convert $plot $plotpng #possible png option???

    echo "finishing plotting $plotpdf"

done < tmp.txt

rm tmp.txt
rm new.cpt
rm $outdir/*.ps


