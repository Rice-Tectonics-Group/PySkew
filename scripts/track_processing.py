#!/usr/bin/env python3
"""
Script which weaves together the many components of picking out data to deskew at a specific magnetic chron for the purpose of automating this process, and provides a number of utilities for visulizing the preprocessing and later analysis of magnetic skewness data.

Flags:
        -h : Displays this help message then exits
    Augmentation Flags:
        -lp   : leave plots open when running -ds, -pbs, -pbsz, or -srp for the sake of reading off data.
    Plot from Deskew Flags:
        -ds   : plot deskew figures shows the current best phase shift and 3 phase shifts a step size away in
               both directions as detailed in the deskew file passed in.
        -pbsz : plot by spreading zone, runs -pbs for each spreading zone seperately and sorts by lattitude.
        -pll  : plot lunes, plots all lunes of given deskew file together.
        -pal  : plot all lunes, plots all lunes of given deskew file seperately.
        -pp   : plot pole, plots all lunes like -pll but also takes a .ellipse file which has lon, lat, angle,
                major, minor axis so it can plot the pole on the same figure.
        -odsk : similar to plot skewness by spreading zone, it plots overlying deskew plots for each spreading zone after the reduction to the pole. 
                It takes two .deskew files and two pole names used to run -rtp flag. 
    Spreading Rate Analysis Flags:
        -srp  : Spreading rate picks, displays your left edge picks of anomaly location for each anomaly in a .srp
              file (see -srf for description of file) next to a synthetic with anomaly location age picks if 
              available. This requires a .deskew file to display the phase shifted profiles followed by a .srp 
              file to display the picks.
        -srr  : Spreading rate results, displays the "results" of your spreading rate analysis by visualizing 
              a .sr file (see -srf for details on creation of .sr file). This displays the mean (or median if -
              median flag is used) of the spreading rate with age with 95% confidence intervals. -title can also 
              be used to force a title on the plot.
        -iso  : This creates a plot of isochrons from a .deskew file and a .srp file (see -srf for details on .srp
              files). This is done by using the left edge picks of anomalies to determine their center and then 
              transform from anomaly data space using the original data files to latlon space where they can be 
              plotted on a map.
    File Generation Flags:
        -rtp  : Reduce to pole, takes a deskew file and a .ellipse file then reduces phase shifts to a pole
              outputing a new deskew file
        -srf  : Create Spreading rate (.sr) file, requires a spreading rate picks file (.srp) consisting of a 
              header with track names and a first column of anomoly names then a list of the picked left edges of
              as many anomalies as possible, as well as a timescale file containing the same anomalies as in 
              the .srp and the column "base" which has the bottom edge of all the anomolies. This then calculates
              the spreading rate mean, median, and std of each anomoly in the two files.
        -srm  : Create spreading rate model file used in generation of the synthetic and calculation of anomalous
              skewness. This generates a file that has the mapping of ages to spreading rates for each spreading 
              zone you are analysing by compiling a series of .sr files generated during your spreading rate 
              analysis. Note that after running this flag you will need to rename the top lines delimenating the 
              different spreading zones so that the names match those in your .deskew file. In addition if you 
              want the code to automatically pull this model you need to place it in this path ../raw_data/
              spreading_rate_model.txt
        -cdmf : Creates a bunch of $DATAFILE.deskewed files which have "lat lon mag" deskewed in it for the phase
               shift given in the input deskew file so it can then be used in gmt scripts for visualization.
        -cor  : Correct intersepts, the automatic intercepts calculated by the geographic preprocessing are often
               off by ~30-50 km and thus should be corrected to make analysis easier and results more accurate.
               This flag does this taking in a deskew file and a correction file.
        -mt   : Create maxtab file, creates a maxtab file from a given deskew file so it can be used with max to
               preform other statistical processing.
        -dt   : Create datatable, deprecated flag creates Lin's matlab datatable file format was origionally used
               to run the matlab code this program depends on and as a check of functionality.
    Geographic Preprocessing Flags:
        -f    : .inp file input, this is REQUIRED for all flags bellow this as it details the basics of what is
               needed to start analysis from the raw data.
        -1    : First stage geographic preprocessing: creates spreading zone files from Barckhausen 2013 and
               Cande 1989 data, finds intersects of all tracks in $CWD/../raw_data/hi_alt and $CWD/../raw_data/ship
               with the chron of interest in the spreading zones just separated out and saves this data to your
               usable_tracks_and_intersects_for_chron#.txt file, separates out the data that is selected as
               intersecting into your data directory, extracts shipmag .latlon file as preprocessing for step 2.
               You MUST run aeromag_prep after this step before step 2. (Note: This process can take some time)
        -2    : Second stage geographic preprocessing: finds all turns in the separated tracks from step 1 using
               RDP and checking angles in the simplified path, cuts tracks and flips according to the heading 
               supplied in your .inp file, find great cicle poles for all spreading zones and save to .gcp files
               for step 3, find which cuts intersect your chron of interest still and saves these to your 
               usable_tracks_and_intersects_for_chron#.txt file, preprocesses shipmag cuts as needed to fit in 
               this framework for analysing their magnetic data. (Note: running -r after this is suggested)
        -3    : Third stage geographic preprocessing: Find azimuth and strike of all spreading zones and 
               intersecting cuts using the intersection point and great circle poles calculated in step 2 then
               save these to .azszs files, plot the geographic data of all cut tracks being kept for analysis so
               they can be evaluated. You MUST run aeromag_prep on all aeromag data after this step.
        -4    : Fourth stage geographic preprocessing: this step puts it all together grabbing from all the 
               previously generated files to create a .deskew file which should be easy to edit and contain all
               relivent information to proceed deskewing the profiles. MUST be run with -a.
        -a    : ages flag, MUST be run with -4. It takes in two ages supplied after the flag which should be the
               min and max age of the anomoly of interest.
        -ptl  : Plot transit locations. This flag plots all intercepts found in the current track data being
               analyzed so that you can visualize the current geographic spread of your dataset.
        -A    : Annotations. This plots annotations on the transit locations plot (only a few spreading zones 
               supported, file may require manual editing to properly move or add annotations). Can only be run
               with -ptl
        -r    : Remove failed data. This does what it says it removes any data that is no longer detailed in the
               usable_tracks_and_intersects_for_chron#.txt file and places it in your failed data directory so it
               can be preserved, but still not cluttering your analysis. It is recommended your run this 
               everytime you backup and edit your usable_tracks_and_intersects_for_chron#.txt file as it will 
               help keep things tidy.

@author: Kevin Gaastra
@date: 5/30/2017
"""

import os,sys,glob
from pyskew.geographic_preprocessing import *
from pyskew.plot_geographic import *
from pyskew.skewness import *
from pyskew.plot_skewness import *

#If run as script it does all processing without plots in one go, with the advantage that this can be run in the background over night or concurently with other processes
if __name__=="__main__":
    #do thing with sys.argv and a input file or something later so you don't have to open the stupid file
    if "-h" in sys.argv: help(__name__); sys.exit()

    kwargs={}

    if "-lp" in sys.argv: kwargs['leave_plots_open']=True

    if "-xl" in sys.argv:
        xli = sys.argv.index('-xl')
        kwargs['xlims'] = sys.argv[xli+1].split(',')

    if "-yl" in sys.argv:
        yli = sys.argv.index('-yl')
        kwargs['ylims'] = sys.argv[yli+1]

    if '-rtp' in sys.argv:
        rtpi = sys.argv.index('-rtp')
        deskew_path = sys.argv[rtpi+1]
        try: pole_lon,pole_lat = float(sys.argv[rtpi+2]),float(sys.argv[rtpi+3])
        except (IndexError,ValueError) as e:
            try: pole_lon,pole_lat = list(map(float,open(sys.argv[rtpi+2],'r').readlines()[0].split()))[0:2]
            except (IOError,OSError,ValueError) as e: print("You must provide a valid lattitude and longitude after the deskew path or a valid ellipse file"); sys.exit()
        reduce_to_pole(deskew_path,pole_lon,pole_lat)
        sys.exit()

    fz_loc_path = None
    if '-fzl' in sys.argv:
        fzli = sys.argv.index('-fzl')
        kwargs['fz_loc_path'] = sys.argv[fzli+1]

    if '-odsk' in sys.argv:
        odski = sys.argv.index('-odsk')
        deskew_path = sys.argv[odski+1]
        deskew_path2 = sys.argv[odski+2]
        if '-pn' in sys.argv:
            pni = sys.argv.index('-pn')
            kwargs['pole_name1'] = sys.argv[pni+1]
            kwargs['pole_name2'] = sys.argv[pni+2]
        overlay_skewness_by_spreading_zone(deskew_path, deskew_path2, **kwargs)
        sys.exit()

    if '-pbsz' in sys.argv:
        pbszi = sys.argv.index('-pbsz')
        deskew_path = sys.argv[pbszi+1]
        ridge_loc_path = None
        if '-rdg' in sys.argv:
            rdgi = sys.argv.index('-rdg')
            kwargs['ridge_loc_path'] = sys.argv[rdgi+1]
        plot_skewness_by_spreading_zone(deskew_path, **kwargs)
        sys.exit()

    if '-ffz' in sys.argv:
        ffzi = sys.argv.index('-ffz')
        deskew_path = sys.argv[ffzi+1]
        fz_directory = os.path.join('raw_data','fracture_zones')
        if '-fd' in sys.argv:
            fdi = sys.argv.index('-fd')
            kwargs['fz_directory'] = sys.argv[fdi+1]
        find_fz_crossings(deskew_path,**kwargs)
        sys.exit()

    if '-ds' in sys.argv:
        dsi = sys.argv.index('-ds')
        deskew_path = sys.argv[dsi+1]
        plot_skewnesses(deskew_path,**kwargs)
        sys.exit()

    if '-cdmf' in sys.argv:
        cdmfi = sys.argv.index('-cdmf')
        deskew_path = sys.argv[cdmfi+1]
        create_deskewed_data_file(deskew_path)
        sys.exit()

    if '-gs' in sys.argv:
        gsi = sys.argv.index('-gs')
        age_min = sys.argv[gsi+1]
        age_max = sys.argv[gsi+2]
        spreading_rate_path = sys.argv[gsi+3]
        age_path = sys.argv[gsi+4]
        generate_synth('both',age_min,age_max,spreading_rate_path,age_path)
        sys.exit()

    if '-pp' in sys.argv:
        ppi = sys.argv.index('-pp')
        deskew_path = sys.argv[ppi+1]
        ellipse_file = sys.argv[ppi+2]
        plot_pole_with_lunes(deskew_path,ellipse_file)
        sys.exit()

    if '-polep' in sys.argv:
        polepi = sys.argv.index('-polep')
        deskew_path = sys.argv[polepi+1]
        ellipse_path = sys.argv[polepi+2]
        if '-dis' in sys.argv:
            disi = sys.argv.index('-dis')
            kwargs['dis'] = float(sys.argv[disi+1])
        plot_pole_perturbations(deskew_path,ellipse_path,**kwargs)
        sys.exit()

    if '-pll' in sys.argv:
        plli = sys.argv.index('-pll')
        deskew_path = sys.argv[plli+1]
        plot_lunes_and_save(deskew_path)
        sys.exit()

    if '-pal' in sys.argv:
        pali = sys.argv.index('-pal')
        deskew_path = sys.argv[pali+1]
        plot_all_lunes_seperate(deskew_path)
        sys.exit()

    if '-mt' in sys.argv:
        mti = sys.argv.index('-mt')
        deskew_path = sys.argv[mti+1]
        anomoly_name = sys.argv[mti+2]
        create_maxtab_file(deskew_path,anomoly_name)
        sys.exit()

    if '-dt' in sys.argv:
        dti = sys.argv.index('-dt')
        print("this dt file type is deprecated and only still exists to allow checking the current code base with the old matlab code, but is not actually required to use any of the python code.")
        deskew_path = sys.argv[dti+1]
        create_matlab_datatable(deskew_path)
        sys.exit()

    if '-cor' in sys.argv:
        cori = sys.argv.index('-cor')
        try:
            deskew_path = sys.argv[cori+1]
            site_cor_path = sys.argv[cori+2]
            correct_site(site_cor_path,deskew_path)
        except IndexError: print("please provide a .deskew file and a site_cor file to correct it with after the -cor flag, aborting")
        sys.exit()

    if '-srp' in sys.argv:
        srpi = sys.argv.index('-srp')
        deskew_path = sys.argv[srpi+1]
        spreading_rate_picks_path = sys.argv[srpi+2]
        plot_spreading_rate_picks(deskew_path,spreading_rate_picks_path,**kwargs)
        sys.exit()

    if '-srr' in sys.argv:
        srri = sys.argv.index('-srr')
        sr_path = sys.argv[srri+1]
        if '-median' in sys.argv: kwargs['median_or_mean']='median'
        if '-title' in sys.argv: kwargs['title']=sys.argv[sys.argv.index('-title')+1]
        plot_spreading_rate_results(sr_path,**kwargs)
        sys.exit()

    if '-srf' in sys.argv:
        srfi = sys.argv.index('-srf')
        spreading_rate_picks_path = sys.argv[srfi+1]
        ages_path = sys.argv[srfi+2]
        create_spreading_rate_file(spreading_rate_picks_path,ages_path)
        sys.exit()

    if '-srm' in sys.argv:
        srmi = sys.argv.index('-srm')
        sr_files = sys.argv[srmi+1:]
        make_spreading_rate_model(sr_files)
        sys.exit()

    if '-iso' in sys.argv:
        isoi = sys.argv.index('-iso')
        deskew_path = sys.argv[isoi+1]
        spreading_rate_picks_path = sys.argv[isoi+2]
        plot_isochron_picks(deskew_path,spreading_rate_picks_path,**kwargs)
        sys.exit()

    if '-sip' in sys.argv:
        sipi = sys.argv.index('-sip')
        deskew_path = sys.argv[sipi+1]
        srp_paths = sys.argv[sipi+2:]
        save_iso_picks(deskew_path,srp_paths)
        sys.exit()

    if '-ru' in sys.argv:
        rui = sys.argv.index('-ru')
        deskew_path = sys.argv[rui+1]
        useble_track_path = sys.argv[rui+2]
        update_useable_tracks_from_deskew(deskew_path,useble_track_path)
        sys.exit()

    if '-f' in sys.argv:
        try:
            fi = sys.argv.index('-f')
            fin = open(sys.argv[fi+1],'r')
            lin = [line.split('\t')[0].rstrip("\n") for line in fin.readlines()]
            chron_to_analyse=lin[0].split(',')
            chron,chron_color = chron_to_analyse
            chron_name = "chron%s"%(str(chron))
            chrons_info = [c.split(',') for c in lin[1].split(';')]
            results_directory = lin[2]
            data_directory = lin[3]
            failed_data_dir = lin[4]
            bounding_lats = list(map(float,lin[5].split(',')))
            bounding_lons = list(map(float,lin[6].split(',')))
            e = float(lin[7])
            tolerance = float(lin[8])
            min_angle = float(lin[9])
            plot = int(lin[10])
            heading = lin[11]
            ship_or_aero = lin[12]
            lon_0=sum(bounding_lons)/2
            lat_0=sum(bounding_lats)/2
        except (IOError, IndexError) as e: print("Problem reading in your input file, please run again with a valid input file"); raise e
    else: print("a valid chron input file must be provided for this script to run anything besides the -ds flag in the commandline please provide one after the -f flag and rerun, aborting"); sys.exit()

    #Run step1, search all files for files that potentally cross chron_to_analyse in the region of interest and seperate that data
    if '-1' in sys.argv:
        spreading_zone_files = seperate_chron_into_spreading_zones(chron_to_analyse)

        #This is very time intensive but only needs to be run once per chron skip if already done for this chron
        if ship_or_aero=='aero':
            tracks = glob.glob('../raw_data/hi_alt/**/*.DAT')
        elif ship_or_aero=='ship':
            tracks = glob.glob('../raw_data/ship/**/*.lp')
        else:
            tracks = glob.glob('../raw_data/ship/**/*.lp')
            tracks += glob.glob('../raw_data/hi_alt/**/*.DAT')

        intersecting_tracks,usable_tracks_path = get_track_intersects(chron_to_analyse, tracks, spreading_zone_files, data_directory=data_directory, bounding_lats=bounding_lats, bounding_lons=bounding_lons, e=e)

        #seperates out the data marked out for analysis above int the data_directory
        tracks, seperated_tracks_path = seperate_data(data_directory, usable_tracks_path)

        #generate shipmag latlon files which are generated by aeromag_prep for aeromag files.
        #This is not strictly necessary, but it makes it a lot easier to fit into the existing framework and
        #I'M LAZY
        shipmag_tracks = list(filter(lambda x: 'ship' in x, tracks))
        shipmag_preprocess(shipmag_tracks)

        aeromag_tracks = list(filter(lambda x: 'aero' in x, tracks))
        aeromag_preprocess(aeromag_tracks)

        #plots all selected plots for analysis on a basic ortho map with selected chrons
        if plot: plot_tracks(chrons_info, results_directory, tracks=tracks, lon_0=lon_0, lat_0=lat_0)

        print("please rerun again with -2 as an arg")

    #find the turns in the tracks set to be analysed, cut and flip those tracks, then calculate best fit great circles for spreading zones to get azimuth and strike for projecting data and decide which cuts intersect which spreading zones closely enough to be analysed, and removes data not going to be used in further analysis
    if '-2' in sys.argv:
        #reads in seperated tracks locations and intercepts if you didn't run the above you must run this
        seperated_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
        seperated_tracks_file = open(seperated_tracks_path,'r')
        seperated_tracks = seperated_tracks_file.readlines()
        seperated_tracks_file.close()
        tracks = [track.split()[0] for track in seperated_tracks]

        #slices tracks using rdp according to tolerance and min_angle, set plot to True to check cuts on plots
        track_cuts = find_track_cuts(tracks, chrons_info, results_directory, tolerance=tolerance, min_angle=min_angle, plot=plot)

        cut_tracks,flipped_data = cut_tracks_and_flip(track_cuts, data_directory, heading=heading)

        #Generates .gcp files or great circle pole files from gmt in order to find strike
        spreading_zone_files = seperate_chron_into_spreading_zones(chron_to_analyse)
        for spreading_zone_file in spreading_zone_files:
            subprocess.check_call('gmt fitcircle %s -L3 > %s'%(spreading_zone_file,spreading_zone_file[:-3]+'gcp'),shell=True)

        #Find and save all of the cuts of the tracks that intersect the chron in the region of interest so that the intersect can be used to generate strike and only the intersecting cuts will be used from now on
        tracks,cut_tracks_path = get_track_intersects(chron_to_analyse, cut_tracks, spreading_zone_files, data_directory=data_directory, bounding_lats=bounding_lats, bounding_lons=bounding_lons, e=e)

        #this time it's litterally just taking the .c# file which is actually a .lp file and copying it so it fits the framework
        shipmag_tracks = list(filter(lambda x: 'ship' in x, tracks))
        shipmag_preprocess(shipmag_tracks)

        aeromag_tracks = list(filter(lambda x: 'aero' in x, tracks))
        aeromag_preprocess(aeromag_tracks)

        print("please rerun with -3 to find strikes and plot cut tracks")

    if '-3' in sys.argv:

        #reads in seperated tracks locations and intercepts if you didn't run the above you must run this
        cut_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
        cut_tracks_file = open(cut_tracks_path,'r')
        cut_tracks = cut_tracks_file.readlines()
        cut_tracks_file.close()
        track_sz_and_inters = [track.split('\t') for track in cut_tracks]

        #Uses spreading zone GCP and track/sz intersect to calculate the azimuth and writes out azsz files to save strike and azimuth as well as intersection point of track and anomoly
        tracks, az_files = generate_az_strike_files(track_sz_and_inters, chron_to_analyse, heading, results_directory, plot=plot)

        #plots all selected plots for analysis on a basic ortho map with selected chrons
        if plot: plot_tracks(chrons_info, results_directory, tracks=tracks, track_dir='all_cut_tracks', lon_0=lon_0, lat_0=lat_0, cuts=True)

        print("please rerun with -4 as an arg to generate a .deskew file and start analysis")

    #pull out East and Vertical componenets of the data into their own sub directories for deskewing analysis
    if '-4' in sys.argv:

        if '-a' in sys.argv:
            ai = sys.argv.index('-a')
            try:
                age_min = float(sys.argv[ai+1])
                age_max = float(sys.argv[ai+2])
            except (IndexError,ValueError) as e:
                print("You must supply both a min and a max age bound after the -a flag both of which must be valid floating point numbers, please rerun with valid ages"); sys.exit()
        else: print("-4 must be run with the -a flag followed by your min age and then your max age bounds for the chron being analysed, please rerun with -a"); sys.exit()

        if age_max<age_min:
            print("Your ages are backwards I'm reversing them, but really? You need to go get a coffee or something.")
            tmp_age = age_max
            age_max = age_min
            age_min = tmp_age

        cut_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
        cut_tracks_file = open(cut_tracks_path,'r')
        cut_tracks = [track.split('\t')[0] for track in cut_tracks_file.readlines()]
        cut_tracks_file.close()

        create_deskew_file(chron_name,results_directory,age_min,age_max,data_directory=data_directory,phase_shift=180,step=60)

    if '-ptl' in sys.argv:
        annotate=False
        if '-A' in sys.argv: annotate=True
        print("centering map on %f %f"%(lon_0,lat_0))
        #plot the chrons on a standard ortho map and plot up the names of some of the more common spreading zones
        plot_transit_locations(chron_to_analyse, chrons_info, results_directory, data_directory=data_directory, lon_0=lon_0,lat_0=lat_0,annotate=annotate)

    if '-r' in sys.argv:
        ri=sys.argv.index('-r')
        try: cut_tracks_path=sys.argv[ri+1]
        except IndexError: cut_tracks_path=''
        #default if no file is provided
        if not os.path.isfile(cut_tracks_path): cut_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
        #moves any created cuts that are not indicated as intersecting the chron of interest
        remove_failed_data(cut_tracks_path, data_directory, failed_data_dir=failed_data_dir)
