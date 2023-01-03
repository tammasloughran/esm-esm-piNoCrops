# Update CABLE vegfrac field and set new tiles to the grid box mean
# Arguments are the old and new dump file and the new vegetation fraction ancillary

import numpy as np, sys, umfile, argparse
from um_fileheaders import *
import netCDF4

parser = argparse.ArgumentParser(description="Update vegetation fractions in dump file")
parser.add_argument('-i', dest='ifile', help='Input UM dump')
parser.add_argument('-o', dest='ofile', help='Output UM dump')
parser.add_argument('-f', dest='fracfile', help='New vegetation fraction (ancillary or netCDF)')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                    default=False, help='verbose output')
args = parser.parse_args()

ntiles = 17

# Get old vegetation fraction from dump file
f = umfile.UMFile(args.ifile)

vlist = []
for k in range(f.fixhd[FH_LookupSize2]):
    ilookup = f.ilookup[k]
    lbegin = ilookup[LBEGIN]
    if lbegin == -99:
        break
    if ilookup[ITEM_CODE]==216:
        a = f.readfld(k)
        vlist.append(a)
assert len(vlist)==ntiles, 'Error - expected %d vegetation classes' % ntiles
old_vegfrac = np.array(vlist)

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
    vegfrac[np.isnan(vegfrac)] = f.missval_r
    d.close()
else:
    # Read the vegetation fraction ancillary
    ffrac = umfile.UMFile(args.fracfile)

    vlist = []
    for k in range(ffrac.fixhd[FH_LookupSize2]):
        ilookup = ffrac.ilookup[k]
        lbegin = ilookup[LBEGIN]
        if lbegin == -99:
            break
        assert ilookup[ITEM_CODE]==216, "Field with unexpected stash code %s" % ilookup[ITEM_CODE]
        a = ffrac.readfld(k)
        vlist.append(a)

    # Create a single array with dimensions [vegtype, lat, lon]
    vegfrac = np.array(vlist)

assert vegfrac.shape[0]==ntiles, 'Error - expected %d vegetation classes' % ntiles

if np.all(old_vegfrac == vegfrac):
    print("Vegetation fields are identical. No output file created")
    sys.exit(0)

# # Check that the masks are identical
# old_mask = (old_vegfrac == f.missval_r)
# new_mask = (vegfrac == f.missval_r)
# if not np.all(old_mask == new_mask):
#     print("Error - land sea masks are different")
#     sys.exit(1)

# Fix all 800 tiled CABLE variables 

g = umfile.UMFile(args.ofile, "w")
g.copyheader(f)

k = 0
while k < f.fixhd[FH_LookupSize2]:
    ilookup = f.ilookup[k]
    lbegin = ilookup[LBEGIN]
    if lbegin == -99:
        break
    if 800 <= ilookup[ITEM_CODE] < 920 and ilookup[ITEM_CODE] not in [883, 884, 885, 887, 888]:
        # These don't seem to be necessary
        # or ilookup[ITEM_CODE] in [230, 234, 236]:
        code = ilookup[ITEM_CODE]
        if args.verbose:
            print("Processing", code)
        a = f.readfld(k)
        vlist = [a]
        # Expect another 16 fields with the same code
        for i in range(1,ntiles):
            ilookup = f.ilookup[k+i]
            if ilookup[ITEM_CODE] != code:
                print("Missing tiled fields with", code, k, i)
                sys.exit(1)
            a = f.readfld(k+i)
            vlist.append(a)
        var = np.array(vlist)
        # Grid box mean
        mean = (var*old_vegfrac).sum(axis=0)
        if var.dtype == np.int:
            # 3 layer snow flag is an integer field
            mean = np.round(mean).astype(np.int)
        # If old fraction was zero and new > 0, set to grid box mean
        var = np.where(np.logical_and(old_vegfrac==0, vegfrac>0), mean, var)
        # Set tiles with new zero fraction to zero
        var[vegfrac==0] = 0.
        if var.dtype != np.int:
            # Scale the new grid box mean to match the old
            newmean = (var*vegfrac).sum(axis=0)
            # Leave points with zero mean alone. This may give a divide by zero warning
            # even though such values are masked by the where.
            scale = np.where(newmean==0.,1.,mean/newmean)
            var *= scale
        for i in range(ntiles):
            g.writefld(var[i],k+i)
        k += ntiles
    elif ilookup[ITEM_CODE]==216:
        # USet the new vegetation fractions
        for i in range(ntiles):
            g.writefld(vegfrac[i],k+i)
        k += ntiles
    else:
        if ilookup[ITEM_CODE] == 30:
            # Save the mask (needed for compression)
            a = f.readfld(k)
            g.writefld(a,k)
            g.mask = a
        else:
            a = f.readfld(k, raw=True)
            g.writefld(a,k, raw=True)
        k += 1

g.close()
