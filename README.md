# PySkew: Tools for Skewness Analysis of Marine Magnetic Anomalies

## How to Install

Download or clone the repository above and install Python3 as well as these libraries:  

- numpy (can be installed with pip, anaconda, or aptitude)
- pandas (can be installed with pip, anaconda, or aptitude)
- matplotlib (can be installed with pip, anaconda, or aptitude)
- mpl\_toolkits.basemap (can be installed with pip, anaconda, or aptitude)
- geographiclib (can be installed with pip, anaconda, or aptitude)
- pmagpy (can be installed with pip or at their github repository [here](https://github.com/PmagPy))
- rdp (can be installed with pip)
- wx (can be installed with pip, anaconda, or aptitude) (optional: required for use of the GUIs) 
- qt (can be installed with pip, anaconda, or aptitude) (optional: if installed this backend will be used on mac for data visuals)

If you are new to Python it is suggested that you install python using the [Anaconda framework](https://www.anaconda.com/download)

Once everything is installed you will need to add the main PySkew directory to your PYTHONPATH environment variable (scripts to automate this process are in development). So if you are running a Mac or Linux based operating system add:
```bash
export PYTHONPATH="$PYTHONPATH:$PATHTOPYSKEW"
```
to your .bashrc file or other init file. If you are running Windows you will need to open a terminal and type:
```dos
setx PYTHONPATH $PATHTOPYSKEW
```
Note that you need to replace $PATHTOPYSKEW with the actual path to the main PySkew directory on your disk not simply copy and paste the above.

## Magnetic Track Data Source

This library is meant to be used with data from the NCEI database found [here](https://maps.ngdc.noaa.gov/viewers/geophysics/) which has ship and aero magnetic trackline data. This data is then passed through a preprocessor and put in raw\_data/ship and raw\_data/hi\_alt respectively. Information on preprocessing can be found [here](https://github.com/Caoimhinmg/PySkew/blob/master/raw\_data/hi_alt/README.md) for aeromag and [here](https://github.com/Caoimhinmg/PySkew/blob/master/raw_data/ship/README.md) for ship.

## Geographic Preprocessing

The first phase of analyzing magnetic trackline data is picking out useful tracks from your raw\_data. The parameters needed for this stage can be found in template.inp which gives an example for analysis of C20r and should be repopulated with your analysis parameters. Picking out useful tracks is then done by looking for intersections with isochrons a number of which are provided in the repository from Cande 1989 and Barckhausen 2013. You can use any isochrons you want from [The Global Seafloor Fabric and Magnetic Lineation Data Base Project](http://www.soest.hawaii.edu/PT/GSFML/) which has an extensive library. Currently taking data from this database and putting it in this library is a little complicated but a processing script to automate this is in progress. This intersection with an isochron is used to calculate site location and azimuth of the magnetic lineation. Trackline data is then segmented using the [RDP line simplification algorithm](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm) and looking for turns in the data exceeding a threshold. This is all then output to a .deskew file for further analysis.  

## Skewness Analysis

Skewness analysis can be preformed by a number of different methods including reducing to iteratively improved pole locations or reduction of individual trackline data. As many of these methods as possible are included in the library and this analysis as well as the above preprocessing are preformed through the use of the track\_processing.py script. For help using this script refer to the in code documentation by running:  

```bash
python track_processing.py -h
```

## Licensing

This code can be freely used, modified, and shared. It is licensed under a 3-clause BSD license. See [LICENSE](https://github.com/Rice-Tectonics-Group/PySkew/blob/master/LICENSE) for details.

## References

- Cande, SC, Labrecque JL, Larson RL, Pittman III WC, Golovchenko X.  1989.  Magnetic lineations of the world's ocean basins. :13-13., Tulsa, OK, United States (USA): Am. Assoc. Pet. Geol., Tulsa, OK

- Barckhausen, U., M. Bagge, and D. S. Wilson (2013), Seafloor spreading anomalies and crustal ages of the Clarion-Clipperton Zone, Mar. Geophys. Res., 1-10.
