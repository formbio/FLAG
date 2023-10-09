#!/usr/bin/env python
"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
    gff_to_genbank.py <GFF annotation file> [<FASTA sequence file> <molecule type>]

 FASTA sequence file: input sequences matching records in GFF. Optional if sequences
   are in the GFF
 molecule type: type of molecule in the GFF file. Defaults to DNA, the most common case.
"""
from __future__ import print_function

import sys
import os

from Bio import SeqIO

from BCBio import GFF


def main(gff_file, fasta_file=None, molecule_type="DNA"):
    out_file = "%s.gb" % os.path.splitext(gff_file)[0]
    if fasta_file:
        fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    else:
        fasta_input = {}
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(_check_gff(_fix_ncbi_id(gff_iter), molecule_type), out_file, "genbank")


def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        yield rec


def _check_gff(gff_iterator, molecule_type):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if "molecule_type" not in rec.annotations:
            rec.annotations["molecule_type"] = molecule_type
        yield _flatten_features(rec)


def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.

    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec


if __name__ == "__main__":
    main(*sys.argv[1:])
