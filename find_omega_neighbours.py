#!/usr/bin/env python

import impact_pdb_database as ipd
from sys import stdout, argv, exit



def usage():
    print('USAGE:   python find_omega_neighbours.py -mass=<float>  (-ccs=<float> | -omega=<float>) [ -db=name_of_pdb_database ] [ -rank=[<int> | all] ] [ -num=[ <int> | all] ] [ -massweight=<float> ] [ -omegaweight=<float> ] [ -h ]')
    print('  Default database    = {:s}'.format(ipd.db_default))
    print('  Default pisa_rank   = {:s}'.format(ipd.rank_default))
    print('  Default massweight  = {:f}'.format(ipd.mw_default))
    print('  Default omegaweight = {:f}'.format(ipd.ow_default))
    
if __name__ == '__main__':

    stdout.write('################ IMPACT PROTEOME SEARCH #################\n')
    stdout.write('#       E G Marklund, M T Degiacomi, C V Robinson       #\n')
    stdout.write('#              A J Baldwin, J L P Benesch               #\n')
    stdout.write('# Collisional cross-sections for structural proteomics  #\n')
    stdout.write('#            Structure (201X), XX(X):XX--XX             #\n')
    stdout.write('#########################################################\n\n')

        
    if len(argv) < 3:
        usage()
        exit(1)

    expect = 'None'

    bDB          = False
    bMass        = False
    bCCS         = False
    bOmega       = False
    bRank        = False
    bNum         = False
    bMassweight  = False
    bOmegaweight = False

    for arg in argv[1:]:
        if arg == '-h':
            usage()
            continue

        
        [flag, value] = arg.split('=')

        if flag == '-db':
            bDB = True
            opt_db = value
            
        elif flag == '-mass':
            bMass = True
            opt_mass = value
            
        elif flag == '-ccs':
            if bOmega:
                stdout.write('-ccs and -omega are mutually exclusive. Choose one.\n')
                exit(2)
                
            bCCS = True
            opt_ccs = float(value)
            
        elif flag == '-omega':
            if bCCS:
                stdout.write('-ccs and -omega are mutually exclusive. Choose one.\n')
                exit(2)

            bOmega = True
            opt_omega = float(value)
            
        elif flag == '-rank':
            bRank = True
            if value.isdigit() or value=='all':
                opt_rank = value
            else:
                print(value)
                print(value.isdigit())
                stdout.write('-rank: Permitted values for the pisa rank are \'all\' or integer numbers\n')
                exit(3)

        elif flag == '-num':
            bNum = True
            if value.isdigit() or value=='all':
                opt_num = value
            else:
                print(value)
                print(value.isdigit())
                stdout.write('-num: Permitted values for the number of neighbours are \'all\' or integer numbers\n')
                exit(4)
                
        elif flag == '-massweight':
            bMassweight = True
            opt_mw = value

        elif flag == '-omegaweight':
            bOmegaweight = True
            opt_ow = value
            
        else:
            stdout.write('Unrecognised option {:s}\n'.format(arg))
            exit(5)

    if not (bMass and (bCCS or bOmega)):
        has_db = 'o'
        if bDB:
            has_db = 'x'

        has_mass = 'o'            
        if bMass:
            has_mass = 'x'
            
        has_ccs = 'o'
        if bCCS:
            has_ccs = 'x'
            
        has_omega = 'o'
        if bOmega:
            has_omega = 'x'
            
        has_rank = 'o'
        if bRank:
            has_rank = 'x'

        has_num = 'o'
        if bNum:
            has_rank = 'x'

        has_mw = 'o'
        if bMassweight:
            has_mw = 'x'
            
        has_ow = 'o'
        if bOmegaweight:
            has_ow = 'x'
        
        stdout.write('Need more arguments: -mass({:s}) (-ccs({:s}) | -omega({:s})) [ -db({:s}) ] [ -rank(:s) ] [ -massweight(:s) ] [ -omegaweight(:s) ]\n'.format(
            has_db, has_mass, has_ccs, has_omega, has_rank, has_num, has_mw, has_ow))
        
        exit(6)

        
    if bMassweight:
        mw = float(opt_mw)
    else:
        mw = ipd.mw_default

    if bOmegaweight:
        ow = float(opt_ow)
    else:
        ow = ipd.ow_default

    if bDB:
        db_name = opt_db
    else:
        db_name = ipd.db_default

        
    stdout.write('##### Probe #####\n')
    
    p = ipd.probe()
    
    p.set_mass(float(opt_mass))
    if bCCS:
        p.set_ccs(float(opt_ccs))
        p.finalise()
    else:
        p.set_omega(float(opt_omega))
        p.set_ccs(ipd.calc_ccs(p.m, p.omega))
        
    p.dump(stdout)

    stdout.write('\nMaking database instance\n')
    db = ipd.pdb_ccs()
    
    stdout.write('\nReading database from file {:s}'.format(db_name))
    
    if bRank:
        rank = opt_rank
        stdout.write(' using pisa rank {:s}\n'.format(rank))
    else:
        stdout.write('\n')
        rank=ipd.rank_default

    if bNum:
        numneighbours = opt_num
        stdout.write(' will print {:s} neighbours\n'.format(numneighbours))
    else:
        stdout.write('\n')
        numneighbours=ipd.num_default
        
    db.read_database(db_name, rank=rank)

    stdout.write('\nFinding neighbours\n\n')
    db.find_neighbours(p, massweight=mw, omegaweight=ow, numneighbours=numneighbours)
    stdout.write('              {:>5s} {:>5s} {:>5s} {:>12s} {:>10s} {:>6s}'.format('Dist.', 'PDB', 'PISA', 'Mass (Da)', 'CCS (A^2)', 'omega'))
    stdout.write('\n')
    db.print_neighbours(stdout)
