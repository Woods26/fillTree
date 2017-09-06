#!/usr/bin/env python

import logging
import os
import sys

from Bio import Phylo
from Bio import SeqIO
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

log = logging.getLogger("phylo-fasta-fudger")


# function definitions
def closest_available_relatives(sp, avail_set, terminals, tree):
    pairs = ((terminals[sp], terminals[c]) for c in avail_set)
    return sorted([(i[1].name, tree.distance(*i)) for i in pairs], key=lambda x: x[1])


def main(fasta_path, output_path, tree_fn):
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
                        for seq_record in SeqIO.parse(fasta_fn[ortho], "fasta")}
    log.debug("Fasta directory loaded. {} fasta files parsed".format(len(fasta)))

    # if output_path does not exist, create it
    if not os.path.isdir(output_path):
        log.debug("Creating output directory {}")
        os.makedirs(output_path)

    # if it does exist, use if empty, or create and use new unique directory
    else:
        output_path = create_unique_dir(output_path)

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


def create_unique_dir(path, limit=99):
    width = len(str(limit))
    original = path.rstrip(os.sep)
    if len(os.listdir(original)) == 0:
        return original  # folder empty, let's use it
    count = 1
    while count < limit:
        try:
            os.mkdir(path)
            log.debug("Using output directory {}".format(path))
            return path
        except OSError as e:
            if e.errno == 17:  # file exists
                path = "{0}_{1:0>{2}}".format(original, count, width)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))


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
    parser.add_argument('-o', '--output_path', help='path to output fasta directory',
                        default="output/")
    parser.add_argument('-t', '--tree_fn', help='path to newick tree file', default="input/tapir_ref_959genes.tre")
    args = parser.parse_args()

    # run main
    try:
        exit(main(args.fasta_path, args.output_path, args.tree_fn))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
