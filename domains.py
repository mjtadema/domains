#!/usr/bin/env python3
from pathlib import Path
from copy import deepcopy

domdef = Path("domains.dat")
assert domdef.exists()

class Domain:
    def __init__(self):
        self.segments = [] # list of segments
        self.exceptions = []

    @property
    def ranges(self):
        for segment in self.segments:
            fr, to = segment
            yield range(fr, to+1)

    def add_segment(self, fr, to):
        self.segments.append((fr, to))

    def __contains__(self, item):
        return any(item in subrange for subrange in self.ranges)

    def __add__(self, l):
        shifted = []
        for (fr, to) in self.segments:
            shifted.append((fr+l, to+l))
        self.segments = shifted
        return self

structure = []
with open(domdef) as fin:
    for line in fin:
        line = line.strip().split()
        n_segments = int(line[0])
        domain = Domain()
        try:
            for _ in range(n_segments):
                line = next(fin)
                fr, to = line.strip().split()
                fr, to = [int(n) for n in (fr, to)]
                domain.add_segment(fr, to)
        except StopIteration:
            pass
        structure.append(domain)

def extrapolate_domains(structure, l, n):
    l = 184
    n = 8
    n_orig_domains = len(structure)
    for i in range(n * n_orig_domains):
        domain_copy = deepcopy(structure[i])
        structure.append(domain_copy + l)
    return structure

structure = extrapolate_domains(structure, 184, 8)

def same_domain(domains: list, fr, to) -> bool:
    for domain in domains:
        if fr in domain and to in domain:
            return True
    return False

def apply_domains(structure):
    atom_to_res = {}

    with open("molecule_0.itp") as fin ,\
    open("molecule_0_doms.itp", 'w') as fout:
        current_header = ""
        for line in fin:
            try:
                if current_header == "atoms":
                    #    1 Q5      1 MET BB     1    1
                    idx, _, resi, *_ = line.strip().split()
                    idx, resi = [int(n) for n in (idx, resi)]
                    atom_to_res[idx] = resi
                elif current_header == "bonds":
                    #    1    3 1 0.350 4000
                    fr, to, bt, l, f = line.strip().split()
                    fr, to = [int(n) for n in (fr, to)]
                    fr, to = [atom_to_res[n] for n in (fr, to)]
                    l, f = [float(n) for n in (l, f)]
                    if f == 1300: # elastic bond force
                        if not same_domain(structure, fr, to):
                            fout.write(";") # comment out offending lines
            except ValueError:
                pass # empty lines probably
            if line.startswith("["):
                current_header = line.strip("[ ]\n")
            fout.write(line)

apply_domains(structure)