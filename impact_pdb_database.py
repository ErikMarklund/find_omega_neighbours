from math import log10, sqrt
from operator import itemgetter

large_dist_const = 100000000000
db_default   = 'pdb_database.dat' 
mw_default   = 1.
ow_default   = 5.
rank_default = '1'
lineparams = [0.457126, 2/3.]

def calc_omega(m, ccs):
    ccs_fit = 10**lineparams[0] * m**lineparams[1]
    omega = ccs/ccs_fit

    return omega

def calc_ccs(m, omega):
    """The inverse of calc_omega"""
    ccs_fit = 10**lineparams[0] * m**lineparams[1]
    ccs = omega * ccs_fit

    return ccs

def calc_omega_mass_distance(m1, omega1, m2, omega2, massweight, omegaweight):
    # Compare the masses on a log-scale.
    return sqrt((massweight * log10(m1/m2))**2 + (omegaweight * (omega1-omega2))**2)


class probe:
    def __init__(self, m=0.0, ccs=0.0):
        self.m = m
        self.ccs = ccs

    def set_mass(self, m):
        self.m = m

    def set_ccs(self, ccs):
        self.ccs = ccs

    def set_omega(self, omega):
        self.omega = omega

    def calc_omega(self):
        try:
            self.omega = calc_omega(self.m, self.ccs)
        except AttributeError:
                print('ERROR: Need either omega or ccs')
                raise

    def finalise(self):
        """Makes sure omega is set"""
        try:
            if self.m <= 0:
                print('ERROR: Probe mass is zero or negative.')
                raise ValueError
                
        except AttributeError:
            print('ERROR: Probe mass unset')
            
        try:
            self.omega += 0
        except AttributeError:
            self.calc_omega()

    def dump(self, fs):
        fs.write('mass  = {:f} Da\n'.format(self.m))
        fs.write('ccs   = {:f} A^2'.format(self.ccs))
        fs.write('\n')
        fs.write('omega = {:f}\n'.format(self.omega))
        
            

class pdb_entry:
    def __init__(self, m=0.0, ccs=0.0, name='', pisa_rank=0, pdistance=large_dist_const, omega=[]):
        self.name = name
        self.m    = m
        self.ccs  = ccs
        if omega==[]:
            self.calc_omega()
        self.pdistance = pdistance
        self.pisa_rank = pisa_rank

    def calc_omega(self):
        self.omega = calc_omega(self.m, self.ccs)

    def set_probedistance(self, p, massweight, omegaweight, mc=5, oc=0.2):
        #p.finalise()
        #if False or (abs(log10(self.m) - log10(p.m))) > mc or (abs(self.omega - p.omega)) > oc:
        #    self.pdistance = large_dist_const
        #else:
        self.pdistance = calc_omega_mass_distance(p.m, p.omega, self.m, self.omega, massweight, omegaweight)
        
    def write(self, fs):
        s = '{:5.3f} {:>5s} {:5d} {:12.2f} {:10.2f} {:6.3f}\n'.format(self.pdistance,
                                                                     self.name,
                                                                     self.pisa_rank,
                                                                     self.m,
                                                                     self.ccs,
                                                                     self.omega)
        fs.write(s)

    def dump(self, fs):
        fs.write('mass      = {:f}\n'.format(float(self.m)))
        fs.write('ccs       = {:f}\n'.format(float(self.ccs)))
        fs.write('omega     = {:f}\n'.format(float(self.omega)))
        fs.write('rank      = {:f}\n'.format(float(self.pisa_rank)))
        
        
class pdb_ccs:
    def __init__(self, entry=[]):
        self.entry = entry
        self.neighbours = []
    
    def add_entry(self, entry):
        self.entry.append(entry)
        
    def calc_omega(self):
        for e in self.entry:
            e.calc_omega()
            
    def read_database(self, fname, rank='1'):

        if rank=='all':
            bSlice = False
        else:
            bSlice = True
            r = int(rank)
        
        with open(fname, 'r') as f:
            try:
                for line in f:
                    stripline = line.strip()
                    if len(stripline) < 1 or stripline[0] in '%#;':
                        continue

                    sline = stripline.split()

                    if len(sline) != 4:
                        print('Badly formatted line')
                        print(line.rstrip())
                        continue

                    name       = sline[0]
                    pisa_rank  = int(sline[1])
                    mass       = float(sline[2])
                    ccs        = float(sline[3])

                    if mass <= 0 or ccs <= 0 or (bSlice and r!=pisa_rank):
                        continue

                    self.add_entry(pdb_entry(m=mass, ccs=ccs, name=name, pisa_rank=pisa_rank))

            except IOError:
                print('Error while reading database.')
                raise

    def get_rank(self, rank):
        db = pdb_ccs()
        for e in self.entry:
            if e.pisa_rank == rank:
                db.add_entry(e)

        return db

    def find_neighbours(self, p, massweight=mw_default, omegaweight=ow_default, mass_cutoff=5, omega_cutoff=0.2):
        self.neighbours = []

        N = []
        
        for e in self.entry:
            e.set_probedistance(p, massweight, omegaweight, mc=mass_cutoff, oc=omega_cutoff)
            N.append((e.name, e.m, e.ccs, e.omega, e.pdistance, e.pisa_rank))
            
        N_sorted = sorted(N, key=itemgetter(4))

        for e in N_sorted[:10]:
            self.neighbours.append(pdb_entry(name=e[0], m=e[1], ccs=e[2], omega=e[3], pdistance=e[4], pisa_rank=e[5]))
        for n in self.neighbours:
            n.calc_omega()
    def print_neighbours(self, fs):
        for (i,n) in enumerate(self.neighbours):
            fs.write('Neighbour {:<2d}: '.format(i))
            n.write(fs)
        
if __name__ == '__main__':
    from sys import stdout

    e = pdb_entry(m=1.00e4, ccs=1.1e3, pisa_rank=1, name='tst')

    e.dump(stdout)
    
    print('Making database instance')
    db = pdb_ccs()
    
    print('Reading database from file')
    db.read_database(ipd.db_default, rank='1')

    print('Making probe')
    p = probe()
    p.set_mass(1e4)
    p.set_ccs(1.1e3)
    p.finalise()
    p.dump(stdout)

    print('Finding neighbours')
    db.find_neighbours(p)

    print('Neighbours:')
    db.print_neighbours(stdout)
    
