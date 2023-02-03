#!/usr/bin/env python
# Update CABLE vegfrac field and set new tiles to the grid box mean.
# Arguments are the old and new dump file and the new vegetation fraction ancillary.
# Note: the variable name vegfrac does not refer to fractions of only vegetated types, but rather
# tile fractions of any type.
# Author: Martin Dix, Tammas Loughran
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


def get_field(restart_file, code):
    """Get field `code` from a UM restart file. Replaces missing values with np.nan.
    """
    var = []
    for k in range(restart_file.fixhd[FH_LookupSize2]):
        ilookup = restart_file.ilookup[k]
        lbegin = ilookup[LBEGIN]
        if lbegin==-99: break
        if ilookup[ITEM_CODE]==code: var.append(restart_file.readfld(k))
    var = np.array(var)
    var[var==restart_file.missval_r] = np.nan
    return var


# Parse arguments.
parser = argparse.ArgumentParser(description="Update vegetation fractions in dump file")
parser.add_argument('-i', dest='ifile', help="Input UM dump")
parser.add_argument('-o', dest='ofile', help="Output UM dump")
parser.add_argument('-f', dest='fracfile', help="New vegetation fraction (ancillary or netCDF)")
parser.add_argument(
        '-n',
        dest='name_fracvar',
        default='field1391',
        help="Name of the fraction variable in the netCDF file",
        )
parser.add_argument(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        default=False,
        help='verbose output',
        )
args = parser.parse_args()

# Get old vegetation fraction from dump file.
f = umfile.UMFile(args.ifile)
old_vegfrac = get_field(f, VEGFRAC_CODE)
old_previous_year = get_field(f, PREV_VEGFRAC_CODE)
assert len(old_vegfrac)==NTILES, 'Error - expected %d vegetation classes' % NTILES
assert len(old_previous_year)==NTILES, 'Error - expected %d vegetation classes' % NTILES

if args.fracfile.endswith(".nc"):
    # Read from a netCDF version of a dump file.
    ncfile = netCDF4.Dataset(args.fracfile)
    v = ncfile.variables[args.name_fracvar]

    # There may be some points outside the valid range.
    v.set_auto_mask(False)
    vegfrac = v[:].squeeze()
    vegfrac = vegfrac.astype(old_vegfrac.dtype)

    # Normalise sums to exactly 1.
    #vegfrac /= vegfrac.sum(axis=0)
    vegfrac[vegfrac==f.missval_r] = np.nan
    ncfile.close()
else:
    # Read the vegetation fraction ancillary
    ffrac = umfile.UMFile(args.fracfile)
    vegfrac = get_field(ffrac, VEGFRAC_CODE)
    ffrac.close()

assert vegfrac.shape[0]==NTILES, 'Error - expected %d vegetation classes' % NTILES

# Check that there are changes in the vegetation fraction for the current and previous year.
if np.all(old_vegfrac==vegfrac)&np.all(old_vegfrac==old_previous_year):
    print("Vegetation fields are identical. No output file created")
    sys.exit(0)

# Check that the masks are identical.
old_mask = (old_vegfrac==f.missval_r)
new_mask = (vegfrac==f.missval_r)
if not np.all(old_mask == new_mask):
    sys.exit("Error - land sea masks are different")

# Update or copy fields to new restart file.
output_file = umfile.UMFile(args.ofile, "w")
output_file.copyheader(f)
k = 0
while k<f.fixhd[FH_LookupSize2]:
    ilookup = f.ilookup[k]
    lbegin = ilookup[LBEGIN]
    if lbegin==-99:
        break

    # Initialize soil and snow temperatures if new tiles were created.
    if ilookup[ITEM_CODE] in [801,802,803,804,805,806,825,826,827]:
        code = ilookup[ITEM_CODE]
        if args.verbose:
            print("Updating", code)
        vlist = [f.readfld(k)]
        for i in range(1, NTILES): # Expect another 16 fields with the same code.
            ilookup = f.ilookup[k+i]
            if ilookup[ITEM_CODE]!=code:
                print("Missing tiled fields with", code, k, i)
                sys.exit(1)
            vlist.append(f.readfld(k+i))
        var = np.array(vlist)

        # Grid box cover fraction weighted mean.
        mean = np.nansum(var*old_vegfrac, axis=0)

        # If old fraction was zero and new>0, set to grid box mean.
        var = np.where(np.logical_and(old_vegfrac==0, vegfrac>0), mean, var)

        # Put missing values back into field and insert into restart file.
        var[old_vegfrac==f.missval_r] = f.missval_r
        var[var==np.nan] = f.missval_r
        for i in range(NTILES):
            output_file.writefld(var[i], k+i)
        k += NTILES

    # If we are resetting the previous year's cover fractions, just use old_vegfrac.
    elif ilookup[ITEM_CODE]==PREV_VEGFRAC_CODE:
        var = old_vegfrac
        var[var==np.nan] = f.missval_r
        for i in range(NTILES):
            output_file.writefld(var[i], k+i)
        k += NTILES

    # Set the new vegetation fractions for the current year.
    elif ilookup[ITEM_CODE]==VEGFRAC_CODE:
        vegfrac[vegfrac==np.nan] = f.missval_r
        for i in range(NTILES):
            output_file.writefld(vegfrac[i], k+i)
        k += NTILES
    else:
        if ilookup[ITEM_CODE]==MASK_CODE:
            # Save the mask (needed for compression).
            mask = f.readfld(k)
            output_file.writefld(mask, k)
            output_file.mask = mask
        else:
            # Copy the remaining restat fields to the new file.
            data = f.readfld(k, raw=True)
            output_file.writefld(data, k, raw=True)
        k += 1

output_file.close()

# List of stashvar codes.
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

