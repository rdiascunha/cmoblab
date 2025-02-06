#!/usr/bin/env python
######################################################################
#
# Python3 script to read Gaussian EET fragment output files and process
#         transition dipole moments, site energies and couplings
#         (dipole and couplings are fixed adopting the orientation of the 
#         vector conneting two atoms in each fragment given as input)
#
# Usage:
#        qmmmpol_analysis_final.py [-h] -i INPUT -o OUTPUT
#
# Arguments:
#  -h, --help            show this help message and exit
#  -f FRAGMENTS, --fra FRAGMENTS
#                        Number of fragments in EET calculation
#  -i INPUT, --in INPUT
#                        Input containing atom indexes for each fragment
#                        (one line per fragment containing the 2 atoms
#                        with a connecting vector alinged with dipole).
#  -o OUTPUT, --out OUTPUT
#                        Output file basename.
#
# Author: Renato D. Cunha, PhD <renatodias@ub.edu>
#
#####################################################################

__author__ = "Renato D. Cunha, PhD"
__email_ = "renatodias@ub.edu"
__version__ = 0.1


# Load modules
#####################################################################

import os
import math
import argparse
import sys
import itertools
import collections

# Global variables
#####################################################################

outfile = ""
infile = ""
nfrag = 0
maxene = 0.0
minene = 0.0
mindip = 0.0
indexes = []
enestfile = []
dipstfile = []
foundfile = []
coupfile = []
foundcoupfile = []

# Functions
#####################################################################

def read_log ( filename ):
    '''Reads LOG file.

    Parameters:
    filename    :  (str) Filename

    """
    Returns:
    None
    '''

    global nfrag,minene,maxene,mindip,indexes,enestfile,dipstfile,foundfile,coupfile,foundcoupfile

    # Store all lines in the LOG file.
    fin = open(filename, "r")
    lines = [line.rstrip('\n') for line in fin]

    # Read lines
    dipst = []
    enest = []
    foundst = []
    coupst = []
    foundcoupst = []
    dip = []
    ene = []
    coords = []
    nfrag1 = 0
    nfrag2 = 0
    nfold = 0
    flagcoup = False
    flagdip  = False
    for line in lines:

        if 'fragment=' and 'resnum' in line.lower():
           line1 = line.replace('=',' ')
           line2 = line1.replace(',',' ')
           # split line into words
           words = line2.split()
           nf = int(words[1])

           if nf != nfold:
               if nfrag1 > 0:
                  coords.append( coordsf )
               nfrag1 += 1
               coordsf = []
               print("Reading coordinates of fragment",nfrag1)

           coordsf.append( [float(words[4]), float(words[5]), float(words[6])] )
           nfold = nf

        if 'Leave Link  101' in line:
            coords.append( coordsf )

        if 'EET: doing calculation for fragment' in line:
            nfrag2 += 1
            dipf = []
            enef = []

        if nfrag2 > 0:
            if 'Ground to excited state transition electric dipole moments (Au)' in line:
               flagdip = True

            if 'Ground to excited state transition velocity dipole moments (Au)' in line:
               flagdip = False

        if nfrag2 > 0 and flagdip:
           if not 'state' in line:
                 # split line into words
                 words = line.split()
                 dipx = float(words[1])
                 dipy = float(words[2])
                 dipz = float(words[3])
                 dipmod = math.sqrt(dipx**2+dipy**2+dipz**2)
                 dipf.append( [dipx, dipy, dipz, dipmod] )

        if nfrag2 > 0 and ' nm  f=' in line:
                 # split line into words
                 words = line.split()
                 enef.append( float(words[4]) )

        if 'Leave Link  914' in line:
            dip.append( dipf )
            ene.append( enef )

    if nfrag1 != nfrag:
           error(' Number of fragments in coordinates does not match number of fragments ')
    if nfrag2 != nfrag:
           error(' Number of EET fragment results does not match number of fragments ')

    for I in range(nfrag):
        print("Coordinates fragment %i" % I)
        print(coords[I])
        print("Dipoles fragment %i" % I)
        print(dip[I])
        print("Energies fragment %i" % I)
        print(ene[I])

    # Identify state of interest and decide dipole orientation
    print("Atomic indexes for fragments")
    print(indexes)

    if len(coords) != nfrag or len(dip) != nfrag or len(ene) != nfrag :
       error(' Number of fragments in coords, dip or ene does not match ')

    state = []
    changesign = []
    for I in range(nfrag):
        for J in range(len(coords[I])):
            if J+1 == indexes[I][0]:
               x1 = coords[I][J][0]
               y1 = coords[I][J][1]
               z1 = coords[I][J][2]

            if J+1 == indexes[I][1]:
               x2 = coords[I][J][0]
               y2 = coords[I][J][1]
               z2 = coords[I][J][2]

        # Dipole expected orientation
        r = []
        r.append( x2 - x1 )
        r.append( y2 - y1 )
        r.append( z2 - z1 )
 
        maxr = 0.0
        maxj = 0
        for J in range(len(r)):
            if abs(r[J]) > abs(maxr):
                maxr = r[J]
                maxj = J

        # Identify state of interest and decide dipole orientation
        print("Atom position vector for fragment %i: %10.3f %10.3f %10.3f" % ( I, r[0], r[1], r[2] ) )
        print("Maximum component: %i %10.3f" % ( maxj, maxr ) )

        if len(dip[I]) != len(ene[I]):
            error(' Different number of states in dipoles and energies ')

        ener = 0.0
        dipmax = -1.0
        state.append(-1)
        for J in range(len(dip[I])):
            if (dip[I][J][3]) > dipmax and ene[I][J] < maxene and (dip[I][J][3]) > mindip and ene[I][J] > minene:
                dipmax = dip[I][J][3]
                state[I] = J

        if state[I] == -1:
            print("Warning: State not found for frag %i in file %s" % (I,filename))
            dipst.append( [0.0, 0.0, 0.0, 0.0] )
            enest.append( 0.0 )
            foundst.append( False )
            changesign.append( False )
        else:
        # Check dipole sign
            if maxr > 0.0 and (dip[I][state[I]][maxj]) > 0.0:
                changesign.append( False )
            elif maxr > 0.0 and (dip[I][state[I]][maxj]) < 0.0:
                changesign.append( True )
            elif maxr < 0.0 and (dip[I][state[I]][maxj]) > 0.0:
                changesign.append( True )
            elif maxr < 0.0 and (dip[I][state[I]][maxj]) < 0.0:
                changesign.append( False )
            else:
                changesign.append( False )

            if changesign[I]:
#                dipst.append( [ -1.0*d for d in dip[I][state[I]] ] )
                dipst.append( [  -1.0*dip[I][state[I]][0], -1.0*dip[I][state[I]][1], -1.0*dip[I][state[I]][2], dip[I][state[I]][3] ] )
            else:
                dipst.append( dip[I][state[I]] )
            enest.append( ene[I][state[I]] )
            foundst.append( True )

    enestfile.append( enest )
    dipstfile.append( dipst )
    foundfile.append( foundst )

    for I in range(nfrag):
        print("Fragment %i State assigned %i" % (I+1,state[I]+1) )

    # Initialize coupling matrices
    Coul = []
    Expl = []
    Totl = []
    Foun = []
    for I in range(nfrag):
        CoulI = []
        ExplI = []
        TotlI = []
        FounI = []
        for J in range(nfrag):
            CoulI.append( 0.0 )
            ExplI.append( 0.0 )
            TotlI.append( 0.0 )
            FounI.append( False )
        Coul.append( CoulI )
        Expl.append( ExplI )
        Totl.append( TotlI )
        Foun.append( FounI )

    # Read electronic coupling values.
    for line in lines:
       if 'Coulomb' in line:
            # split line into words
            words = line.split()
            fi = int(words[3]) 
            si = int(words[5].replace(')','')) 
            fj = int(words[8]) 
            sj = int(words[10].replace(')','')) 
            coup = float(words[11]) 
            if state[fi-1] == si - 1 and state[fj-1] == sj - 1:
               Coul[fi-1][fj-1] = coup
               Coul[fj-1][fi-1] = coup
               Foun[fi-1][fj-1] = True
               Foun[fj-1][fi-1] = True

       if 'Explicit' in line:
            # split line into words
            words = line.split()
            fi = int(words[4]) 
            si = int(words[6].replace(')','')) 
            fj = int(words[9]) 
            sj = int(words[11].replace(')','')) 
            coup = float(words[12]) 
            if state[fi-1] == si - 1 and state[fj-1] == sj - 1:
               Expl[fi-1][fj-1] = coup
               Expl[fj-1][fi-1] = coup
               Foun[fi-1][fj-1] = True
               Foun[fj-1][fi-1] = True

       if 'TOTAL COUPLING' in line:
            # split line into words
            words = line.split()
            fi = int(words[4]) 
            si = int(words[6].replace(')','')) 
            fj = int(words[9]) 
            sj = int(words[11].replace(')','')) 
            coup = float(words[12]) 
            if state[fi-1] == si - 1 and state[fj-1] == sj - 1:
               Totl[fi-1][fj-1] = coup
               Totl[fj-1][fi-1] = coup
               Foun[fi-1][fj-1] = True
               Foun[fj-1][fi-1] = True

    for I in range(nfrag):
        for J in range(nfrag):
            if (changesign[I] and not changesign[J] ) or ( changesign[J] and not changesign[I]):
                Coul[I][J] = -Coul[I][J]
                Coul[J][I] = -Coul[J][I]
                Expl[I][J] = -Expl[I][J]
                Expl[J][I] = -Expl[J][I]
                Totl[I][J] = -Totl[I][J]
                Totl[J][I] = -Totl[J][I]

    print("Coul matrix")
    for I in range(nfrag):
        print(Coul[I])

    print("Expl matrix")
    for I in range(nfrag):
        print(Expl[I])

    print("Totl matrix")
    for I in range(nfrag):
        print(Totl[I])

    print("Foun matrix")
    for I in range(nfrag):
        print(Foun[I])

    coupfile.append( [Coul, Expl, Totl] )
    foundcoupfile.append( Foun )
    fin.close()

#
# Standard error:
#

def error(string,where=None):
  # Determine calling function
  if where is None: where = sys._getframe().f_back.f_code.co_name
  header = ' ERROR in %s '  % where
  print("\n %s \n%s\n" % (header.center(50,'-'),string))
  sys.exit()

def parser():
    """ Parses the arguments from the command line and updates the global variables.

    ( ) --> None

    Usage:
    parser()
    """
    global outfile,infile,nfrag,maxene,mindip,minene

    parser = argparse.ArgumentParser(description="ALL.")

    parser.add_argument('-f','--fra',type=int,required=True,help="Number of EET fragments (required).")
    parser.add_argument('-e','--ema',type=float,default=3.5,help="Maximum energy for state of interest (default 3.5).")
    parser.add_argument('-n','--emi',type=float,default=2.0,help="Minimum energy for state of interest (default 2.0).")
    parser.add_argument('-d','--dmi',type=float,default=4.0,help="Minimum dipole for state of interest (default 4.0).")
    parser.add_argument('-i', '--inp',type=str,default="indexes.dat",help="Atom indexes file (default indexes.dat).")
    parser.add_argument('-o', '--out',type=str,default="results",help="Output file basname (default results).")

    args = parser.parse_args()

    # Update global variables
    outfile   = args.out             # Output file basname
    infile    = args.inp             # Input file
    nfrag     = args.fra             # Number of fragments
    mindip    = args.dmi             # Minimum dipole for state
    maxene    = args.ema             # Maximum energy for state
    minene    = args.emi             # Minimum energy for state

# Main
#####################################################################

if __name__ == "__main__":

    # Parser the command line arguments 
    parser()

    print( "Reading Gaussian log files in directory" )
    print( "=======================================" )

    # Store all lines in the file.
    lines = [line.rstrip('\n') for line in open(infile, "r")]

    # Read lines
    nf = 0
    for line in lines:

        # split line into words
        words = line.split()

        # Read and store information
        indexes.append( [int(words[0]), int(words[1])] )
        print("Atom indexes for fragment %i are: %s %s" % ( nf+1, words[0], words[1] ))
        nf += 1

    if nf != nfrag:
       error(' Number of lines in indexes file does not match number of fragments ')

    # Open files for writing average quantities
    foutfrag = []
    fouten = open( outfile+"_ene.dat", "w" )
    foutdi = open( outfile+"_dip.dat", "w" )
    foutcoul = open( outfile+"_coul.dat", "w" )
    foutexpl = open( outfile+"_expl.dat", "w" )
    fouttotl = open( outfile+"_totl.dat", "w" )
    foutscre = open( outfile+"_scre.dat", "w" )
    for I in range(nfrag):
        fname = open( "%s_frag%i.dat" % (outfile,I+1), "w" )
        foutfrag.append( fname )

    # Read the LOG files
    directory = os.getcwd()
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".log"):
            print("Reading file %s\n" % filename)
            read_log(filename)
            continue
        else:
            continue

    # Compute average properties and write all outputs
    # Initialize lists
    Enav   = []
    Dipav  = []
    Nav    = []
    Coulav = []
    Explav = []
    Totlav = []
    Screav = []
    Ncoupav = []
    for I in range(nfrag):
        Enav.append( 0.0)
        Dipav.append( [ 0.0, 0.0, 0.0, 0.0] )
        Nav.append( 0)
        CoulI = []
        ExplI = []
        TotlI = []
        ScreI = []
        NcoupI = []
        for J in range(nfrag):
            CoulI.append( 0.0 )
            ExplI.append( 0.0 )
            TotlI.append( 0.0 )
            ScreI.append( 0.0 )
            NcoupI.append( 0 )
        Coulav.append( CoulI )
        Explav.append( ExplI )
        Totlav.append( TotlI )
        Screav.append( ScreI )
        Ncoupav.append( NcoupI )

    print("Ncoupav initial matrix")
    for I in range(nfrag):
        print(Ncoupav[I])

# Check fragment matrices
    for F in range(len(enestfile)):

        print("Coul matrix file %i" % (F) )
        for I in range(nfrag):
            print(coupfile[F][0][I])

        print("Expl matrix file %i" % (F) )
        for I in range(nfrag):
            print(coupfile[F][1][I])

        print("Totl matrix file %i" % (F) )
        for I in range(nfrag):
            print(coupfile[F][2][I])

        print("Foun matrix file %i" % (F) )
        for I in range(nfrag):
            print(foundcoupfile[F][I])

    for I in range(nfrag):
        foutfrag[I].write(" Energy DipX DipY DipZ Dip\n")

    for F in range(len(enestfile)):
        for I in range(nfrag):
            foutfrag[I].write("%9.3f %9.3f %9.3f %9.3f %9.3f\n" % (enestfile[F][I],dipstfile[F][I][0],dipstfile[F][I][1],dipstfile[F][I][2],dipstfile[F][I][3]) )
            if foundfile[F][I]:
                Enav[I] += enestfile[F][I]
                for K in range(4):
                    Dipav[I][K] += dipstfile[F][I][K]
                Nav[I] += 1
            for J in range(nfrag):
                if foundcoupfile[F][I][J]:
                    Coulav[I][J] += coupfile[F][0][I][J]
                    Explav[I][J] += coupfile[F][1][I][J]
                    Totlav[I][J] += coupfile[F][2][I][J]
                    Screav[I][J] += coupfile[F][2][I][J]/coupfile[F][0][I][J]
                    Ncoupav[I][J] += 1
        print("Ncoupav matrix after file %i" % (F) )
        for I in range(nfrag):
            print(Ncoupav[I])

    for I in range(nfrag):
        if Nav[I] == 0:
            error(' State could never be assigned for a given fragment ')
        Enav[I] = Enav[I]/Nav[I]
        Dipavnew = [ x/Nav[I] for x in Dipav[I]]
        Dipav[I] = Dipavnew
        print("Fragment %i identified in %i frames: %i percent" % (I+1,Nav[I],(Nav[I]/len(enestfile))*100) )
        for J in range(nfrag):
            if I != J:
                if Ncoupav[I][J] == 0 :
                    error(' States could never be assigned for a given pair of fragments ')
                Coulav[I][J] = Coulav[I][J]/Ncoupav[I][J]
                Explav[I][J] = Explav[I][J]/Ncoupav[I][J]
                Totlav[I][J] = Totlav[I][J]/Ncoupav[I][J]
                Screav[I][J] = Screav[I][J]/Ncoupav[I][J]

    for I in range(nfrag):
        fouten.write("%i %9.3f\n" % (I+1,Enav[I]) )
        foutdi.write("%9.3f %9.3f %9.3f %9.3f\n" % (Dipav[I][0], Dipav[I][1], Dipav[I][2], Dipav[I][3]) )
        for J in range(nfrag):
            foutcoul.write("%9.1f" % (Coulav[I][J]) )
            foutexpl.write("%9.1f" % (Explav[I][J]) )
            fouttotl.write("%9.1f" % (Totlav[I][J]) )
            foutscre.write("%9.2f" % (Screav[I][J]) )

        foutcoul.write("\n")
        foutexpl.write("\n")
        fouttotl.write("\n")
        foutscre.write("\n")

    # Close output files
    fouten.close()
    foutdi.close()
    foutcoul.close()
    foutexpl.close()
    fouttotl.close()
    foutscre.close()
    for I in range(nfrag):
        foutfrag[I].close()

