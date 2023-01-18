#!/usr/bin/env python
# Update CABLE vegfrac field and set new tiles to the grid box mean.
# Arguments are the old and new dump file and the new vegetation fraction ancillary.
# Note: the variable name vegfrac does not refer to fractions of only vegetated types, but rather
# tile fractions of any type.
import argparse
import sys

import netCDF4
import numpy as np
import umfile
from um_fileheaders import *

NTILES = 17
VEGFRAC_CODE = 216
PREV_VEGFRAC_CODE = 835
MASK_CODE = 30

# Parse arguments.
parser = argparse.ArgumentParser(description="Update vegetation fractions in dump file")
parser.add_argument('-i', dest='ifile', help='Input UM dump')
parser.add_argument('-o', dest='ofile', help='Output UM dump')
parser.add_argument('-f', dest='fracfile', help='New vegetation fraction (ancillary or netCDF)')
parser.add_argument(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        default=False,
        help='verbose output',
        )
args = parser.parse_args()

# Get old vegetation fraction from dump file
f = umfile.UMFile(args.ifile)
old_vegfrac = []
old_previous_year = []
for k in range(f.fixhd[FH_LookupSize2]):
    ilookup = f.ilookup[k]
    lbegin = ilookup[LBEGIN]
    if lbegin==-99:
        break
    if ilookup[ITEM_CODE]==VEGFRAC_CODE:
        old_vegfrac.append(f.readfld(k))
    if ilookup[ITEM_CODE]==PREV_VEGFRAC_CODE:
        old_previous_year.append(f.readfld(k))
assert len(old_vegfrac)==NTILES, 'Error - expected %d vegetation classes' % NTILES
old_vegfrac = np.array(old_vegfrac)
old_previous_year = np.array(old_previous_year)

if args.fracfile.endswith(".nc"):
    # Read from a netCDF version of a dump file
    d = netCDF4.Dataset(args.fracfile)
    v = d.variables['field1391']
    # There may be some points outside the valid range
    v.set_auto_mask(False)
    vegfrac = v[0]
    vegfrac = vegfrac.astype(old_vegfrac.dtype)
    # Normalise sums to exactly 1
    vegfrac /= vegfrac.sum(axis=0)
    vegfrac[old_vegfrac==f.missval_r] = f.missval_r
    d.close()
else:
    # Read the vegetation fraction ancillary
    ffrac = umfile.UMFile(args.fracfile)
    vegfrac = []
    for k in range(ffrac.fixhd[FH_LookupSize2]):
        ilookup = ffrac.ilookup[k]
        lbegin = ilookup[LBEGIN]
        if lbegin==-99:
            break
        assert ilookup[ITEM_CODE]==VEGFRAC_CODE, "Field with unexpected stash code %s" % ilookup[ITEM_CODE]
        vegfrac.append(ffrac.readfld(k))
    # Create a single array with dimensions [vegtype,lat,lon]
    vegfrac = np.array(vegfrac)

assert vegfrac.shape[0]==NTILES, 'Error - expected %d vegetation classes' % NTILES

if np.all(old_vegfrac==vegfrac)&np.all(old_vegfrac==old_previous_year):
    print("Vegetation fields are identical. No output file created")
    sys.exit(0)

# # Check that the masks are identical
# old_mask = (old_vegfrac == f.missval_r)
# new_mask = (vegfrac == f.missval_r)
# if not np.all(old_mask == new_mask):
#     print("Error - land sea masks are different")
#     sys.exit(1)

# Fix all 800 tiled CABLE variables
output_file = umfile.UMFile(args.ofile, "w")
output_file.copyheader(f)
k = 0
while k < f.fixhd[FH_LookupSize2]:
    ilookup = f.ilookup[k]
    lbegin = ilookup[LBEGIN]
    if lbegin==-99:
        break
    if 800<=ilookup[ITEM_CODE]<920 and ilookup[ITEM_CODE] not in [883,884,885,887,888]:
        # Initialize soil quantities and CNP pools if new tiles were created.
        code = ilookup[ITEM_CODE]
        if args.verbose:
            print("Processing", code)
        vlist = [f.readfld(k)]
        # Expect another 16 fields with the same code
        for i in range(1, NTILES):
            ilookup = f.ilookup[k+i]
            if ilookup[ITEM_CODE]!=code:
                print("Missing tiled fields with", code, k, i)
                sys.exit(1)
            vlist.append(f.readfld(k+i))
        var = np.array(vlist)
        # Grid box cover fraction weighted mean.
        mean = (var*old_vegfrac).sum(axis=0)
        if var.dtype==np.int:
            # 3 layer snow flag is an integer field
            mean = np.round(mean).astype(np.int)
        # If old fraction was zero and new > 0, set to grid box mean
        var = np.where(np.logical_and(old_vegfrac==0, vegfrac>0), mean, var)
        # Set tiles with new zero fraction to zero
        var[vegfrac==0] = 0.0
        # Put missing values back into field
        var[old_vegfrac==f.missval_r] = f.missval_r
        if ilookup[ITEM_CODE]==PREV_VEGFRAC_CODE:
            # If we are resetting the previous year's cover fractions, just use old-vegfrac.
            var = old_vegfrac
        for i in range(NTILES):
            output_file.writefld(var[i], k+i)
        k += NTILES
    elif ilookup[ITEM_CODE]==VEGFRAC_CODE:
        # Set the new vegetation fractions
        for i in range(NTILES):
            output_file.writefld(vegfrac[i], k+i)
        k += NTILES
    else:
        if ilookup[ITEM_CODE]==MASK_CODE:
            # Save the mask (needed for compression)
            mask = f.readfld(k)
            output_file.writefld(mask, k)
            output_file.mask = mask
        else:
            data = f.readfld(k, raw=True)
            output_file.writefld(data, k, raw=True)
        k += 1

output_file.close()

# List of stashvar codes
#s03i801:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 1" ;
#s03i802:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 2" ;
#s03i803:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 3" ;
#s03i804:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 4" ;
#s03i805:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 5" ;
#s03i806:long_name = "CABLE SOIL TEMPERATURE ON TILES LAYER 6" ;
#s03i807:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 1" ;
#s03i808:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 2" ;
#s03i809:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 3" ;
#s03i810:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 4" ;
#s03i811:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 5" ;
#s03i812:long_name = "CABLE SOIL MOISTURE ON TILES LAYER 6" ;
#s03i813:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 1" ;
#s03i814:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 2" ;
#s03i815:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 3" ;
#s03i816:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 4" ;
#s03i817:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 5" ;
#s03i818:long_name = "CABLE FROZEN SOIL MOIST FRAC ON TILES LAYER 6" ;
#s03i819:long_name = "CABLE SNOW DEPTH ON TILES LAYER 1" ;
#s03i820:long_name = "CABLE SNOW DEPTH ON TILES LAYER 2" ;
#s03i821:long_name = "CABLE SNOW DEPTH ON TILES LAYER 3" ;
#s03i822:long_name = "CABLE SNOW MASS ON TILES LAYER 1" ;
#s03i823:long_name = "CABLE SNOW MASS ON TILES LAYER 2" ;
#s03i824:long_name = "CABLE SNOW MASS ON TILES LAYER 3" ;
#s03i825:long_name = "CABLE SNOW TEMPERATURE ON TILES LAYER 1" ;
#s03i826:long_name = "CABLE SNOW TEMPERATURE ON TILES LAYER 2" ;
#s03i827:long_name = "CABLE SNOW TEMPERATURE ON TILES LAYER 3" ;
#s03i828:long_name = "CABLE SNOW DENSITY ON TILES LAYER 1" ;
#s03i829:long_name = "CABLE SNOW DENSITY ON TILES LAYER 2" ;
#s03i830:long_name = "CABLE SNOW DENSITY ON TILES LAYER 3" ;
#s03i831:long_name = "CABLE MEAN SNOW DENSITY ON TILES" ;
#s03i832:long_name = "CABLE SNOW AGE ON TILES" ;
#s03i833:long_name = "CABLE SNOW FLAG ON TILES" ;
#s03i835:long_name = "PREVIOUS YEAR SURF FRACTIONS (TILES)" ;
#s03i851:long_name = "CARBON POOL LABILE ON TILES" ;
#s03i852:long_name = "CARBON POOL PLANT - LEAF ON TILES" ;
#s03i853:long_name = "CARBON POOL PLANT - WOOD ON TILES" ;
#s03i854:long_name = "CARBON POOL PLANT - ROOT ON TILES" ;
#s03i855:long_name = "CARBON POOL LITTER - METB ON TILES" ;
#s03i856:long_name = "CARBON POOL LITTER - STR ON TILES" ;
#s03i857:long_name = "CARBON POOL LITTER - CWD ON TILES" ;
#s03i858:long_name = "CARBON POOL SOIL - MIC ON TILES" ;
#s03i859:long_name = "CARBON POOL SOIL - SLOW ON TILES" ;
#s03i860:long_name = "CARBON POOL SOIL - PASS ON TILES" ;
#s03i861:long_name = "NITROGEN POOL PLANT - LEAF ON TILES" ;
#s03i862:long_name = "NITROGEN POOL PLANT - WOOD ON TILES" ;
#s03i863:long_name = "NITROGEN POOL PLANT - ROOT ON TILES" ;
#s03i864:long_name = "NITROGEN POOL LITTER - METB ON TILES" ;
#s03i865:long_name = "NITROGEN POOL LITTER - STR ON TILES" ;
#s03i866:long_name = "NITROGEN POOL LITTER - CWD ON TILES" ;
#s03i867:long_name = "NITROGEN POOL SOIL - MIC ON TILES" ;
#s03i868:long_name = "NITROGEN POOL SOIL - SLOW ON TILES" ;
#s03i869:long_name = "NITROGEN POOL SOIL - PASS ON TILES" ;
#s03i870:long_name = "NITROGEN POOL SOIL MINIMUM (TILES)" ;
#s03i871:long_name = "PHOSPHORUS POOL PLANT - LEAF (TILES)" ;
#s03i872:long_name = "PHOSPHORUS POOL PLANT - WOOD (TILES)" ;
#s03i873:long_name = "PHOSPHORUS POOL PLANT- ROOT (TILES)" ;
#s03i874:long_name = "PHOSPHORUS POOL LITTER - METB (TILES)" ;
#s03i875:long_name = "PHOSPHORUS POOL LITTER - STR (TILES)" ;
#s03i876:long_name = "PHOSPHORUS POOL LITTER - CWD (TILES)" ;
#s03i877:long_name = "PHOSPHORUS POOL SOIL - MIC (TILES)" ;
#s03i878:long_name = "PHOSPHORUS POOL SOIL - SLOW (TILES)" ;
#s03i879:long_name = "PHOSPHORUS POOL SOIL - PASS (TILES)" ;
#s03i880:long_name = "PHOSPHORUS POOL SOIL LABILE (TILES)" ;
#s03i881:long_name = "PHOSPHORUS POOL SOIL SORB ON TILES" ;
#s03i882:long_name = "PHOSPHORUS POOL SOIL OCC ON TILES" ;
#s03i884:long_name = "NITROGEN DEPOSITION" ;
#s03i885:long_name = "NITROGEN FIXATION" ;
#s03i893:long_name = "LEAF AREA INDEX (CASA-CNP GLAI)" ;
#s03i895:long_name = "WOOD FLUX CARBON (CASA-CNP)" ;
#s03i896:long_name = "WOOD FLUX NITROGEN (CASA-CNP)" ;
#s03i897:long_name = "WOOD FLUX PHOSPHOR (CASA-CNP)" ;
#s03i898:long_name = "WOOD HARVEST CARBON1(CASA-CNP)" ;
#s03i899:long_name = "WOOD HARVEST CARBON2(CASA-CNP)" ;
#s03i900:long_name = "WOOD HARVEST CARBON3(CASA-CNP)" ;
#s03i901:long_name = "WOOD HARVEST NITROG1(CASA-CNP)" ;
#s03i902:long_name = "WOOD HARVEST NITROG2(CASA-CNP)" ;
#s03i903:long_name = "WOOD HARVEST NITROG3(CASA-CNP)" ;
#s03i904:long_name = "WOOD HARVEST PHOSPH1(CASA-CNP)" ;
#s03i905:long_name = "WOOD HARVEST PHOSPH2(CASA-CNP)" ;
#s03i906:long_name = "WOOD HARVEST PHOSPH3(CASA-CNP)" ;
#s03i907:long_name = "WOOD RESPIRA CARBON1(CASA-CNP)" ;
#s03i908:long_name = "WOOD RESPIRA CARBON2(CASA-CNP)" ;
#s03i909:long_name = "WOOD RESPIRA CARBON3(CASA-CNP)" ;
#s03i910:long_name = "WOOD RESPIRA NITROG1(CASA-CNP)" ;
#s03i911:long_name = "WOOD RESPIRA NITROG2(CASA-CNP)" ;
#s03i912:long_name = "WOOD RESPIRA NITROG2(CASA-CNP)" ;
#s03i913:long_name = "WOOD RESPIRA PHOSPH1(CASA-CNP)" ;
#s03i914:long_name = "WOOD RESPIRA PHOSPH2(CASA-CNP)" ;
#s03i915:long_name = "WOOD RESPIRA PHOSPH3(CASA-CNP)" ;
#s03i916:long_name = "THIN RATIO FOR FOREST (CASA-CNP)" ;
#s03i917:long_name = "NITROGEN NET RELEASE (CASA-CNP)" ;
#s03i918:long_name = "NITROGEN LEACHING (CASA-CNP)" ;
#s03i919:long_name = "NITROGEN UPTAKE (CASA-CNP)" ;
#s03i920:long_name = "NITROGEN LOSS (CASA-CNP)" ;

