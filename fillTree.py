#!/usr/bin/env python

import logging
import os
import shutil
import sys

from Bio import Phylo
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


# logging configuration
class OneLineExceptionFormatter(logging.Formatter):
    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        result = super().format(record)
        if record.exc_text:
            result = result.replace(r"\n", "")
        return result


handler = logging.StreamHandler()
formatter = OneLineExceptionFormatter(logging.BASIC_FORMAT)
handler.setFormatter(formatter)
root = logging.getLogger()
root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
root.addHandler(handler)

log = logging.getLogger("fillTree")


# function definitions
def closest_available_relatives(sp, avail_set, terminals, tree):
    pairs = ((terminals[sp], terminals[c]) for c in avail_set)
    return sorted([(i[1].name, tree.distance(*i)) for i in pairs], key=lambda x: x[1])


def main(fasta_path, output_path, tree_fn, remove_old=False):
    log.debug("Parsing newick tree at {}".format(tree_fn))
    tree = Phylo.read(tree_fn, 'newick')
    tree.ladderize()  # Flip branches so deeper clades are displayed at top

    sp_list = [cl.name for cl in tree.get_terminals()]
    log.info("species_list={}".format(str(sp_list)))
    sp_set = set(sp_list)
    terminals = {l.name: l for l in tree.get_terminals()}

    log.debug("Parsing fasta directory {}".format(fasta_path))
    # create handles for all .fasta files in fasta directory
    fasta_fn = {name.split('.')[0]: os.path.join(fasta_path, name) for name in
                os.listdir(fasta_path) if
                ((".fasta" in name) and (".fai" not in name))}

    # read and parse fasta files for each species
    fasta = {}
    for ortho in fasta_fn.keys():
        fasta[ortho] = {seq_record.id: seq_record
                        for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                      "fasta", alphabet=IUPAC.ambiguous_dna)}
    log.debug("Fasta directory loaded. {} fasta files parsed".format(len(fasta)))

    # if output_path does not exist, create it
    if not os.path.isdir(output_path):
        log.debug("Creating output directory {}")
        os.makedirs(output_path)

    # if it exists and remove_old is set, remove and recreate it
    elif remove_old:
        log.debug("Removing and recreating output directory {}")
        shutil.rmtree(output_path)
        os.makedirs(output_path)

    for ortho in fasta.keys():
        avail_set = set(fasta[ortho].keys())
        missing_set = sp_set - avail_set
        file_name = ortho + ".padded.fasta"
        file_path = os.path.join(output_path, file_name)
        if os.path.isfile(file_path):
            log.debug("File {} exists, overwriting".format(file_name))
        log.info("file={}".format(file_name))
        with open(file_path, "wt") as f:
            for sp in sp_list:
                if sp in missing_set:
                    relatives = closest_available_relatives(sp, avail_set, terminals, tree)
                    sub_sp, distance = relatives[0]
                    sub_seq = fasta[ortho][sub_sp].seq
                    description = "subbed-missing-with_{} distance_{}".format(sub_sp, distance)
                    seq_rec = SeqRecord(sub_seq, id=sp, description=description)
                    f.write(seq_rec.format("fasta"))
                    second_sp, second_distance = relatives[1]
                    log.info(
                        "missing={} - substituting={} - distance={} - next_nearest={} - next_distance={}".format(sp,
                                                                                                                 sub_sp,
                                                                                                                 distance,
                                                                                                                 second_sp,
                                                                                                                 second_distance))
                else:
                    f.write(fasta[ortho][sp].format("fasta"))


if __name__ == "__main__":
    import argparse


    # argparse configuration
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)


    parser = MyParser()
    parser.add_argument('-f', '--fasta_path', help='path to input fasta directory', default="input/orthoExon_fasta/")
    parser.add_argument('-o', '--output_path', help='path to output fasta directory (be careful, will delete contents',
                        default="output/")
    parser.add_argument('-t', '--tree_fn', help='path to newick tree file', default="input/tapir_ref_959genes.tre")
    parser.add_argument('-r', '--remove_old', help='If flag is set, OUTPUT_PATH will be deleted and recreated.',
                        action='store_true')
    args = parser.parse_args()

    # run main
    try:
        exit(main(args.fasta_path, args.output_path, args.tree_fn, args.remove_old))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
