#!/bin/bash
#ADD SITE LOCATIONS TO MAP TO DISPLAY CURRENT ISOCHRONS
#In file defs
infile=$1
plot="`dirname $infile`/magwiggle.ps"
plotpdf="`dirname $infile`/magwiggle.pdf"
grdfile=/Users/alice/Documents/GitHub/PySkew/ETOPO1_Bed_g_gmt4.grd

#Map boundary defs
y1=-30
y2=25
x1=230
x2=270

#Flags
region=-R$x1/$x2/$y1/$y2
proj=-Jm0.2c
cpt=-Cnew.cpt
cptfile=new.cpt
zscale=-Z600

#Script start

gmt makecpt -Cpolar -T-5000/-2000/150 > $cptfile

 grdimage $grdfile $cpt $proj $region -K > $plot

 pscoast $region $proj -W1p,black -Ggrey -B5 -K -O >> ${plot}

 psscale -Dn1.2/0+w10/1 -O -K $region $proj $cpt >> ${plot}

while read i; do
    if [ "${i:0:9}" != "comp_name" ]; then
        echo $i | awk '{if (substr($1,1,1) != "#") printf("%s\t%s/%s.deskewed\t%s\n", $11, $9, $1, $8-180)}' >> tmp.txt
        echo $i | awk '{if (substr($1,1,1) != "#") printf("%s,%s\n", $7, $6)}' >> tmp_latlon.txt
    fi
done < $infile

while read t f azi; do
    if [ $t == "aero" ]; then
        color="black"
    else
        color="darkslategrey"
    fi
     psxy tmp_latlon.txt $region $proj -O -K -Glimegreen -Sc.1 -t75 >> ${plot}
     pswiggle $f $region $proj $zscale -O -K -I$azi -G-$color -Wthin,$color -Tthin,$color,- -t50 >> ${plot}
done < tmp.txt

rm tmp.txt
rm tmp_latlon.txt
rm new.cpt

echo "saving $plotpdf"

ps2pdf $plot $plotpdf

rm $plot

#xdg-open $plotpdf #hack remove when done debugging
open $plotpdf #hack remove when done debugging


