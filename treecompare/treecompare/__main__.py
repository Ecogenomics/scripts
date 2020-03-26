###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import argparse
import os
import sys

import luigi

from treecompare import __version__
from treecompare.pipeline.compare_trees import CompareTreesPipeline
from treecompare.pipeline.f_measure_tree import GenerateFMeasureTree


def print_help():
    lines = ['',
             f'...::: treecompare v{__version__} :::...'.center(80),
             '',
             '\tfmeasure',
             '\t\tCreate a dendrogram showing ranks which are not fully monophyletic.',
             '',
             '\tcompare',
             '\t\tCompare the topological/taxonomic differences between two or more trees.',
             '']
    print('\n'.join(lines))


def assert_has_args(argparser, lst_args):
    missing_args = list()
    for arg in lst_args:
        if not hasattr(argparser, arg) or getattr(argparser, arg) is None:
            missing_args.append(arg)

    if len(missing_args) > 0:
        print('you are missing the following arguments:')
        print(', '.join(missing_args))
        sys.exit(1)


def assert_has_files(paths):
    missing_paths = list()
    for path in paths:
        if not os.path.isfile(path):
            missing_paths.append(path)

    if len(missing_paths) > 0:
        print('the following paths do not exist:')
        print('\n'.join(missing_paths))


def main(args=None):
    # Create the argument parsers.
    parser = argparse.ArgumentParser(description='comparison of phylogenetic '
                                                 'tree topology and taxonomic '
                                                 'differences between trees')
    subparsers = parser.add_subparsers(help='--', dest='subparser_name')

    # Create the f-measure tree subparser.
    subparser_fmt = subparsers.add_parser('fmeasure',
                                          help='create a dendrogram showing ranks which are not fully monophyletic')
    subparser_fmt.add_argument('--trees', type=str, help='input tsv containing trees compare')
    subparser_fmt.add_argument('--taxonomy', type=str, help='file indicating taxonomy of each genome')
    subparser_fmt.add_argument('--out_dir', type=str, help='path to store the results in')
    subparser_fmt.add_argument('--cpus', type=int, default=1, help='number of cpus to use')

    # Create the tree comparison subparser
    subparser_cmp = subparsers.add_parser('compare',
                                          help='compare the topological/taxonomic differences between two or more trees')

    subparser_cmp.add_argument('--trees', type=str, help='input tsv containing trees compare')
    subparser_cmp.add_argument('--taxonomy', type=str, help='file indicating taxonomy of each genome')
    subparser_cmp.add_argument('--out_dir', type=str, help='path to store the results in')
    subparser_cmp.add_argument('--common', type=str, default='global',
                               help='use {global, local} common taxa when comparing differences')
    subparser_cmp.add_argument('--cpus', type=int, default=1, help='number of cpus to use')

    # Verify that a subparser was selected
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'treecompare v{__version__}')
        sys.exit(0)
    else:
        print(f'treecompare v{__version__}')
        args = parser.parse_args(args)

        # Determine which pipeline to execute
        if args.subparser_name == 'fmeasure':
            assert_has_args(args, ('trees', 'taxonomy', 'out_dir', 'cpus'))
            assert_has_files((args.trees, args.taxonomy))
            luigi.build([GenerateFMeasureTree(path_trees=args.trees,
                                              dir_out=args.out_dir,
                                              path_taxonomy=args.taxonomy)],
                        local_scheduler=True, workers=max(1, args.cpus))

        elif args.subparser_name == 'compare':
            assert_has_args(args, ('trees', 'taxonomy', 'out_dir', 'cpus', 'common'))
            assert_has_files((args.trees, args.taxonomy))
            luigi.build([CompareTreesPipeline(path_trees=args.trees,
                                              dir_out=args.out_dir,
                                              path_taxonomy=args.taxonomy,
                                              common=args.common)],
                        local_scheduler=True, workers=max(1, args.cpus))

        else:
            print_help()
            sys.exit(0)

    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()
