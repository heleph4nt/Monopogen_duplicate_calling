"""
The MIT License
...
"""

# This file provides utilities for constructing, writing, parsing,
# filtering, and comparing VCF files.
#
# It is largely inherited from MonoVar and reused by Monopogen
# to standardize VCF output and post-processing.

import argparse
import os
import sys
import gzip
from subprocess import Popen, PIPE, check_call
import re
import time
from glob import glob


# Template for the VCF header metadata section
VCF_meta_template = """##fileformat=VCFv4.1
##fileDate={_t.tm_year}-{_t.tm_mon}-{_t.tm_mday}
##source=MonoVar
{_d.FILTER_META}{_d.INFO_META}{_d.FORMAT_META}{_d.REF_META}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{_d.FILES_META}
"""

# Template for a single VCF record line
VCF_record_template = "{_r.CHROM}\t{_r.POS}\t{_r.ID}\t{_r.REF}\t{_r.ALT}\t{_r.QUAL}\t{_r.FILTER}\t{_r.INFO}\t{_r.FORMAT}\t{_r.PASSCODE}\n"


class VCFRecord:
    # Lightweight container for INFO fields of a VCF record
    def __init__(self):
        self.info = {}


class VCFDocument:
    # This class manages writing a complete VCF file:
    # header metadata, reference information, and variant records.

    def __init__(self, outf):
        self.time = time.ctime()
        self.info_fields = []
        self.filter_fields = []
        self.format_fields = []
        self.outf = outf

    def populate_fields(self, bam_id_list):
        # Define standard FILTER, FORMAT, and INFO fields
        # used by MonoVar / Monopogen

        self.filter_fields.append(('LowQual', 'Low quality'))

        self.format_fields.append(
            ('AD', '.', 'Integer', 'Allelic depths for the ref and alt alleles'))
        self.format_fields.append(
            ('DP', '1', 'Integer', 'Approximate read depth'))
        self.format_fields.append(
            ('GQ', '1', 'Integer', 'Genotype Quality'))
        self.format_fields.append(
            ('GT', '1', 'String', 'Genotype'))
        self.format_fields.append(
            ('PL', 'G', 'Integer', 'Phred-scaled genotype likelihoods'))

        self.info_fields.append(
            ('AC', 'A', 'Integer', 'Allele count in genotypes'))
        self.info_fields.append(
            ('AF', 'A', 'Float', 'Allele frequency'))
        self.info_fields.append(
            ('AN', '1', 'Integer', 'Total number of alleles'))
        self.info_fields.append(
            ('BaseQRankSum', '1', 'Float', 'Base quality rank-sum test'))
        self.info_fields.append(
            ('DP', '1', 'Integer', 'Read depth'))
        self.info_fields.append(
            ('QD', '1', 'Float', 'Quality by depth'))
        self.info_fields.append(
            ('SOR', '1', 'Float', 'Strand bias odds ratio'))
        self.info_fields.append(
            ('MPR', '1', 'Float', 'Log odds ratio of non-ref observation'))
        self.info_fields.append(
            ('PSARR', '1', 'Float', 'Alt-to-ref read ratio per sample'))

        # Sample (cell) IDs
        self.files_list = bam_id_list

    def populate_reference(self, ref_file):
        # Register reference genome used for calling
        self.ref_file = ref_file

    def print_header(self):
        # Write the full VCF header to file

        self.FILTER_META = ''.join(
            '##FILTER=<ID=%s,Description="%s">\n' % _ for _ in self.filter_fields)
        self.FORMAT_META = ''.join(
            '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.format_fields)
        self.INFO_META = ''.join(
            '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.info_fields)

        self.FILES_META = '\t'.join(self.files_list)
        self.REF_META = '##reference=file:{0}\n'.format(self.ref_file)

        self.outf.write(VCF_meta_template.format(_d=self, _t=time.localtime()))

    def print_record(self, record):
        # Write a single variant record with INFO fields expanded
        record.INFO = ';'.join(
            "%s=%s" % (_[0], str(record.info[_[0]]))
            for _ in self.info_fields if _[0] in record.info
        )
        self.outf.write(VCF_record_template.format(_r=record))

    def print_my_record(self, record):
        # Write a preformatted VCF record directly
        self.outf.write(VCF_record_template.format(_r=record))

    def close(self):
        self.outf.close()


class VRecord:
    # Represents a single VCF variant record

    def __init__(self, chrm, pos):
        self.CHROM = chrm
        self.POS = pos

    def get6fields(self, ref, alt, id, qual, filter_i, info):
        # Populate the core VCF fields (ID, REF, ALT, QUAL, FILTER, INFO)
        self.ID = id
        self.REF = ref
        self.ALT = alt
        self.QUAL = str(qual)
        self.FILTER = filter_i
        self.INFO = info

    def get_passcode(self, barcode):
        # Store per-sample genotype string(s)
        self.PASSCODE = barcode

    def format_vcf(self, samples):
        # Define FORMAT field and ordering of per-sample columns
        self.FORMAT = "\t".join(samples)


class VCF:
    # Utility class for reading and querying existing VCF files

    def __init__(self, fn):
        self.fn = fn
        self.header = ''
        self.chrmline = ''

    def open(self):
        # Open plain or gzipped VCF
        if self.fn.endswith(".gz"):
            self.fh = gzip.open(self.fn, "r")
        else:
            self.fh = open(self.fn, "r")

    def close(self):
        self.fh.close()

    def read_header(self):
        # Read header and extract sample list
        while True:
            line = self.fh.readline()
            if not line:
                break
            if line.startswith("#CHROM"):
                self.samplelist = line.strip().split()[9:]
                self.chrmline = '\t'.join(line.split()[:9])
                break
            else:
                self.header += line

    def format_header(self):
        # Return reconstructed header line
        return self.header + self.chrmline + '\t' + '\t'.join(self.samplelist)

    def read1(self):
        # Generator over VCF records
        while True:
            line = self.fh.readline()
            if not line:
                break
            if line[0] == '#':
                continue
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.id = pair[2]
            r.data = dict(
                zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]])
            )
            yield r

    def fetch_region(self, chrm, beg, end):
        # Query a genomic region using tabix
        for line in Popen(
            ["tabix", self.fn, "%s:%d-%d" % (chrm, beg, end)],
            stdout=PIPE
        ).stdout:
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.data = dict(
                zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]])
            )
            yield r


def main_compare(args):
    # Compare genotype calls between two VCF files
    # across shared samples and variant IDs

    vcf1 = VCF(sys.argv[1])
    vcf2 = VCF(sys.argv[2])

    vcf1.open()
    vcf2.open()

    vcf1.read_header()
    vcf2.read_header()

    intersample = set(vcf2.samplelist) & set(vcf1.samplelist)
    print(len(intersample), "overlapping samples")

    cs2g2 = {}
    for r in vcf2.read1():
        for s in intersample:
            cs2g2[(r.id, s)] = r.data[s]

    cs2g1 = {}
    for r in vcf1.read1():
        for s in intersample:
            cs2g1[(r.id, s)] = r.data[s]

    overlap = set(cs2g1.keys()) & set(cs2g2.keys())

    # Confusion matrix of genotype concordance
    vars = ["./.", "0/0", "0/1", "1/1"]
    for var1 in vars:
        sys.stdout.write(var1)
        for var2 in vars:
            sys.stdout.write('\t%d' % len(
                [_ for _ in overlap if cs2g1[_] == var1 and cs2g2[_] == var2]
            ))
        sys.stdout.write('\n')


def main_filter(args):
    # Filter VCF records by variant ID list

    vcf = VCF(args.v)
    vcf.open()
    vcf.read_header()

    cids = set()
    with open(args.cid) as fh:
        for line in fh:
            pair = line.split('\t')
            cids.add(pair[args.cidcol - 1])

    print(vcf.format_header())
    for r in vcf.read1():
        if r.id in cids:
            print(vcf.format1(r))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='vcf tool')
    subparsers = parser.add_subparsers()

    psr_compare = subparsers.add_parser("compare", help="compare vcfs")
    psr_compare.add_argument('-v1', help='VCF file 1')
    psr_compare.add_argument('-v2', help='VCF file 2')
    psr_compare.set_defaults(func=main_compare)

    psr_filter = subparsers.add_parser("filter", help="filter vcf")
    psr_filter.add_argument('-v', help='VCF file')
    psr_filter.add_argument('--cid', help='call id list')
    psr_filter.add_argument('--cidcol', type=int, default=1)
    psr_filter.set_defaults(func=main_filter)

    args = parser.parse_args()
    args.func(args)
