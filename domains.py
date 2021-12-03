#!/usr/bin/env python3
"""
Script to make domains in a course grained elastic network.
Domains are defined by a domdef.dat file.

Written by Matthijs J. Tadema, MSc (2020)
"""
from pathlib import Path
from copy import deepcopy
import argparse
import logging
logging.basicConfig(level=logging.INFO)


class Domain:
    def __init__(self):
        self.segments = [] # list of segments
        self.exceptions = []
        self.bonds = []

    @property
    def ranges(self):
        for segment in self.segments:
            fr, to = segment
            yield range(fr, to+1)

    def add_segment(self, fr, to):
        self.segments.append((fr, to))

    def add_exception(self, fr, to):
        self.exceptions.append({fr, to})

    def add_bonds(self , fr, to, l, f, bt=1):
        self.bonds.append((fr, to, bt, l, f))

    def __contains__(self, item):
        return any(item in subrange for subrange in self.ranges)

    def __add__(self, l):
        "move all the atom indices by l"
        shifted = []
        for (fr, to) in self.segments:
            shifted.append((fr+l, to+l))
        self.segments = shifted
        excepted = []
        for (fr, to) in self.exceptions:
            excepted.append({fr+l, to+l})
        self.exceptions = excepted
        bonds = []
        for (fr, to, *rest) in self.bonds:
            bonds.append((fr+l, to+l, *rest))
        self.bonds = bonds
        return self

    def shift(self, l: int, start: int, end: int, nshift: int):
        """
        Shift part of the definition by nshift
        :param l: int, length of the protein
        :param start: int
        :param end: int
        :param nshift: int
        :return: None
        """
        shifted = []
        modfix = lambda n: ((n - 1) % l) + 1
        shiftrange = range(start, end+1)
        for (fr, to) in self.segments:
            if fr in shiftrange:
                fr += nshift
            if to in shiftrange:
                to += nshift
            shifted.append((fr, to))
        self.segments = shifted
        excepted = []
        for (fr, to) in self.exceptions:
            if modfix(fr) in shiftrange:
                fr += nshift
            if modfix(to) in shiftrange:
                to += nshift
            # This was the crucial step to fix making bonds between chains
            fr += (fr // l) * nshift
            to += (to // l) * nshift
            excepted.append({fr, to})
        self.exceptions = excepted
        bonds = []
        for (fr, to, *rest) in self.bonds:
            if modfix(fr) in shiftrange:
                fr += nshift
            if modfix(to) in shiftrange:
                to += nshift
            # This was the crucial step to fix making bonds between chains
            fr += (fr // l) * nshift
            to += (to // l) * nshift
            bonds.append((fr, to, *rest))
        self.bonds = bonds

    def wrap(self, l, n):
        shifted = []
        total_l = l*n
        # This is because residues are 1 indexed,
        # but mod goes to zero...
        modfix = lambda n: ((n-1) % total_l) + 1
        for (fr, to) in self.segments:
            shifted.append((modfix(fr), modfix(to)))
        self.segments = shifted
        excepted = []
        for (fr, to) in self.exceptions:
            excepted.append({modfix(fr), modfix(to)})
        self.exceptions = excepted
        bonds = []
        for (fr, to, *rest) in self.bonds:
            bonds.append((modfix(fr), modfix(to), *rest))
        self.bonds = bonds


def parse_range(line):
    fr, to = line.strip().split()
    fr, to = [int(n) for n in (fr, to)]
    return fr, to


def parse_domains(domdef: Path) -> list:
    """
    Parse a domdef file and return list of domains.
    :param domdef: Path
    :return: list
    """
    domdef = Path(domdef)
    assert domdef.exists()
    structure = []
    with open(domdef) as fin:
        for line in fin:
            if line.startswith("#"): continue
            line = line.strip().split()
            n_segments = int(line[0])
            try:
                n_exceptions = int(line[1])
            except IndexError:
                n_exceptions = 0
            try:
                n_new_bonds = int(line[2])
            except IndexError:
                n_new_bonds = 0
            logging.debug(f"{n_segments} segments, {n_exceptions}, {n_new_bonds} new bonds.")
            domain = Domain()
            try:
                for _ in range(n_segments):
                    line = next(fin)
                    fr, to = parse_range(line)
                    logging.debug(f"Segment between {fr} {to}")
                    domain.add_segment(fr, to)
                for _ in range(n_exceptions):
                    line = next(fin)
                    domain.add_exception(*parse_range(line))
                for _ in range(n_new_bonds):
                    line = next(fin)
                    fr, to, l, f = line.strip().split()
                    fr, to = [int(n) for n in (fr, to)]
                    l, f = [float(n) for n in (l, f)]
                    domain.add_bonds(fr, to, l, f)

            except StopIteration:
                pass
            structure.append(domain)
    return structure


def extrapolate_domains(structure: list, *, length: int, nchains: int, start=None, end=None, nshift=0, **kwargs) -> list:
    """
    Extrapolate given domains to [n] domains, give length [l]
    :param structure: list
    :param length: int
    :param nchains: int
    :param start: int
    :param end: int
    :param nshift: int
    :return: list, new structure
    """
    if not any(obj is None for obj in (start, end, nshift)):
        for domain in structure:
            # Shift segments and bonds
            # if part of the structure was shifted
            domain.shift(length, start, end, nshift)

    n_orig_domains = len(structure)
    for i in range((nchains - 1) * n_orig_domains):
        # Extrapolate domains to several chains
        domain_copy = deepcopy(structure[i])
        structure.append(domain_copy + length)
    for domain in structure:
        # Finally wrap the indices around if they exceed mod total length
        domain.wrap(length, nchains)
    return structure


def same_domain(domains: list, fr, to) -> bool:
    for domain in domains:
        if fr in domain and to in domain:
            return True
    return False


def is_exception(domains, fr, to):
    for domain in domains:
        if {fr, to} in domain.exceptions:
            return True
    return False


def apply_domains(structure, *, itp_in: Path, itp_out: Path, fc: float, **kwargs):
    """
    Apply the given domain definitions to an itp file
    (i.e. trim the elastic network to be only in domains)
    :param structure: list
    :param itp_in: Path
    :param itp_out: Path
    :return:
    """
    atom_to_res = {}
    res_to_bb = {}

    counter = 0

    logging.info(f"Writing to {itp_out.name}")
    with open(itp_in) as fin ,\
    open(itp_out, 'w') as fout:
        current_header = ""
        for line in fin:
            try:
                if current_header == "atoms":
                    #    1 Q5      1 MET BB     1    1
                    idx, _, resi, resn, name, *_ = line.strip().split()
                    idx, resi = [int(n) for n in (idx, resi)]
                    atom_to_res[idx] = resi
                    if name == "BB":
                        res_to_bb[resi] = idx
                elif current_header == "bonds":
                    #    1    3 1 0.350 4000
                    fr, to, bt, l, f = line.strip().split()
                    fr, to = [int(n) for n in (fr, to)]
                    fr, to = [atom_to_res[n] for n in (fr, to)]
                    l, f = [float(n) for n in (l, f)]
                    if not is_exception(structure, fr, to):
                        if f == fc: # elastic bond force
                            if not same_domain(structure, fr, to):
                                counter += 1
                                logging.debug(f"Cutting {fr} {to}")
                                fout.write(";") # comment out offending lines
                    else:
                        logging.debug(f"Found exception! {fr} {to}")

            except ValueError:
                pass # empty lines probably
            if line.startswith("["):
                current_header = line.strip("[ ]\n")
                if current_header == "constraints":
                    # Start writing new bonds
                    for domain in structure:
                        for bond in domain.bonds:
                            fr, to, bt, l, f = bond
                            fr, to = [res_to_bb[n] for n in (fr, to)]
                            fout.write(f"{fr} {to} {bt} {l} {f}\n")
                    fout.write("\n")
            fout.write(line)
    if counter == 0:
        logging.warning(f"only cut {0} bonds, check domain definitions or fc (was {fc}).")


def main(*, domdef, length: int, nchains: int, **kwargs):
    logging.info("Parsing domains")
    structure = parse_domains(domdef)
    if not any(obj is None for obj in (length, nchains)):
        logging.info("Extrapolating domains")
        structure = extrapolate_domains(structure, length=length, nchains=nchains, **kwargs)
    logging.info("Applying domains")
    apply_domains(structure, **kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="Be loud and noisy.")
    parser.add_argument("--force-constant", "-fc", default=500, type=float, dest="fc",
                        help="Force constant to detect elastic bonds to cut.")
    files = parser.add_argument_group('files')
    files.add_argument("--domdef", type=Path, default="domdef.dat",
                        help="Domain definition (default: domdef.dat)")
    files.add_argument("--itp-in", "-i", default="molecule_0.itp", type=Path,
                        help="Input topology (default: molecule_0.itp)")
    files.add_argument("--itp-out", "-o", default="molecule_0_doms.itp", type=Path,
                        help="Output topology (default: molecule_0_doms.itp)")
    extra = parser.add_argument_group('extrapolation')
    extra.add_argument("--nchains", "-n", type=int, help="Number of (identical) chains in the topology")
    extra.add_argument("--length", "-l", type=int, help="Length of a single chain")
    shift = parser.add_argument_group('shift')
    shift.add_argument("--start", type=int)
    shift.add_argument("--end", type=int)
    shift.add_argument("--nshift", type=int, default=0)

    args = parser.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Logging level debug")
    assert args.domdef.exists()

    main(**vars(args))