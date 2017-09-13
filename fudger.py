#!/usr/bin/env python
# coding=utf-8
"""Create fasta files with the complete set of species found in a reference tree by padding with closest relatives.

If importing as a module, see the function 'main' for usage documentation.
If running as a script, run with the -h flag for usage documentation.
"""

import logging
import os
import random
import sys
import argparse

from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# logging configuration
# noinspection PyMissingOrEmptyDocstring
class OneLineExceptionFormatter(logging.Formatter):
    def formatException(self, exc_info):
        result = super().formatException(exc_info)
        return repr(result)

    def format(self, record):
        """Format error message to fit on a single line.
        """
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
def closest_available_relatives(sp, avail_sp_set, terminals_dict, tree):
    """Return a list of (species, distance) tuples from avail_sp_set sorted by distance to target species.

    :param sp: Species of interest, typically a missing species for which a substitute must be found
    :param avail_sp_set: Subset of species in terminals_dict which are not missing (e.g. from a fasta file)
    :param terminals_dict: A dictionary of terminal nodes from tree keyed by name (i.e. species)
    :param tree: A phylogenetic tree object from the Bio.Phylo package
    :return: A list of (species, distance) tuples from avail_sp_set sorted by distance to target species sp
    """
    pairs = ((terminals_dict[sp], terminals_dict[c]) for c in avail_sp_set)
    return sorted([(i[1].name, tree.distance(*i)) for i in pairs], key=lambda x: (x[1], random.random()))


def parse_fasta_file(file_name):
    """Return a dictionary of records from the given fasta file indexed by their record id

    :param file_name: The path to the target fasta file
    :return: a dictionary of fasta records indexed by their record id
    """
    return {seq_record.id: seq_record
            for seq_record in SeqIO.parse(file_name, "fasta")}


def parse_fasta_filenames_from_directory(fasta_path):
    """ Return a dictionary of paths to fasta files in the directory fasta_path indexed by the section of file name
    preceding the first '.' (i.e. the ortholog)

    :param fasta_path: Path to a directory containing fasta files
    :return: A dictionary of paths to fasta files indexed by ortholog
    """
    return {name.split('.')[0]: os.path.join(fasta_path, name) for name in
            os.listdir(fasta_path) if
            ((".fasta" in name) and (".fai" not in name))}


def create_unique_dir(path, limit=99):
    """Return a path to an empty directory. Either the dir at path, or a dir of the form 'path + _01'

    :param path: The initial path to use
    :param limit: The maximum number of directory variations this function will attempt to create.
    :return: A path to an empty directory.
    """
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
        except OSError as path_error:
            if path_error.errno == 17:  # file exists
                path = "{0}_{1:0>{2}}".format(original, count, width)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))


def write_fasta_file_with_substitutetions(ortho, fasta, sp_list, output_path, terminals_dict, tree):
    """Write the fasta file for ortho substituting the sequence of the closest relative for any missing species.

    :param ortho: Target ortholog
    :param fasta: Dict of fasta records for ortho indexed by record id (i.e. species)
    :param sp_list: List of all species represented in tree. Order determines order of fasta entries in output file
    :param output_path: Path to new output fasta file
    :param terminals_dict: A dictionary of terminal nodes from tree keyed by name (i.e. species)
    :param tree: A phylogenetic tree object from the Bio.Phylo package
    """
    avail_sp_set = set(fasta.keys())
    missing_sp_set = set(sp_list) - avail_sp_set
    out_file_name = ortho + ".padded.fasta"
    out_file_path = os.path.join(output_path, out_file_name)
    if os.path.isfile(out_file_path):
        log.debug("File {} exists, overwriting".format(out_file_name))
    log.info("file={}".format(out_file_name))
    with open(out_file_path, "wt") as f:
        for sp in sp_list:
            if sp not in missing_sp_set:
                f.write(fasta[sp].format("fasta"))
            else:
                relatives_list = closest_available_relatives(sp, avail_sp_set, terminals_dict, tree)
                substitute_sp, distance = relatives_list[0]
                substitute_seq = fasta[substitute_sp].seq
                description = "subbed-missing-with_{} distance_{}".format(substitute_sp, distance)
                seq_rec = SeqRecord(substitute_seq, id=sp, description=description)
                f.write(seq_rec.format("fasta"))
                second_sp, second_distance = relatives_list[1]
                # noinspection PyPep8
                log.info(
                    "missing={} - substituting={} - distance={} - next_nearest={} - next_distance={}".format(sp,
                                                                                                             substitute_sp,
                                                                                                             distance,
                                                                                                             second_sp,
                                                                                                             second_distance))


def main(fasta_path, output_path, tree_fn, seed=None, sort_tree=False):
    # set seed for random function optionally using user supplied seed for repeatability
    """Parse the newick tree and directory of fasta files and output a new directory of fasta
    files with closest relatives substituted for any missing species.

    :param fasta_path: Input fasta directory path
    :param output_path: Output fasta directory path
    :param tree_fn: Path to the newick tree file
    :param seed: Optional seed value used to make random function repeatable
    :param sort_tree: Optional bool value determining whether to sort the tree. Determines species order of output files
    """
    if seed is None:
        seed = random.randrange(sys.maxsize)
    random.seed(seed)
    log.info("random_seed={}".format(seed))

    # parse tree from file
    log.debug("Parsing newick tree at {}".format(tree_fn))
    tree = Phylo.read(tree_fn, 'newick')

    # optionally sort the tree
    if sort_tree:
        tree.ladderize()  # Flip branches so deeper clades are displayed at top

    # parse names from tree
    sp_list = [cl.name for cl in tree.get_terminals()]
    log.info("species_list={}".format(str(sp_list)))

    # create terminals dictionary indexed by name
    terminals_dict = {leaf.name: leaf for leaf in tree.get_terminals()}

    # create handles for all .fasta files in fasta directory
    log.debug("Parsing fasta directory {}".format(fasta_path))
    fasta_fn_dict = parse_fasta_filenames_from_directory(fasta_path)

    # read and parse fasta files for each species
    fasta_dict = {ortho: parse_fasta_file(path) for ortho, path in fasta_fn_dict.items()}
    log.debug("Fasta directory loaded. {} fasta files parsed".format(len(fasta_dict)))

    # if output_path does not exist, create it
    if not os.path.isdir(output_path):
        log.debug("Creating output directory {}")
        os.makedirs(output_path)
    else:  # If it exists and is empty, use it. If it exists and is not empty, create and use new unique directory
        output_path = create_unique_dir(output_path)

    for ortho, fasta in fasta_dict.items():
        write_fasta_file_with_substitutetions(ortho, fasta, sp_list, output_path, terminals_dict, tree)


if __name__ == "__main__":

    # argparse configuration
    # noinspection PyMissingOrEmptyDocstring
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
    parser.add_argument('-r', '--repeatable_seed', help='optionally set seed to make random repeatable', default=None)
    parser.add_argument('-s', '--sort_tree', help='optionally sort the tree so deeper clades are displayed at top',
                        action='store_true')
    args = parser.parse_args()

    # run main
    try:
        exit(main(args.fasta_path, args.output_path, args.tree_fn, seed=args.repeatable_seed, sort_tree=args.sort_tree))
    except Exception as e:
        log.exception("Exception in main(): {}".format(e))
        exit(1)
