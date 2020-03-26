import os
import subprocess
from collections import defaultdict

import dendropy
import luigi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from treecompare.external.genetreetk import GeneTreeTk
from treecompare.io import InputTreesFile, NewickTreeFile, TaxonomyFile
from treecompare.tools import maybe_make_dirs, convert_tax_to_ranks

DIR_INT = 'intermediate_results'
DIR_DECORATE = os.path.join(DIR_INT, 'decorate')

DIR_OUT = 'results'


class DecorateTree(luigi.Task):
    """Decorates a rooted tree using PhyloRank.

            return DecorateTree(
            input_tree='/home/uqamussi_local/github/treecompare/example/trees/00.ar122_r89_reference_rooted.tree',
            taxonomy_file=self.path_taxonomy,
            output_tree=os.path.join(self.dir_out, DIR_DECORATE, 'foo', 'foo.tree'))

    """
    tree_id = luigi.Parameter()
    input_tree = luigi.Parameter()
    taxonomy_file = luigi.Parameter()
    output_tree = luigi.Parameter()

    def requires(self):
        return {'input_tree': NewickTreeFile(path=self.input_tree),
                'taxonomy_file': TaxonomyFile(path=self.taxonomy_file)}

    def run(self):
        maybe_make_dirs(self.output_tree)
        args = ['phylorank', 'decorate', self.requires()['input_tree'].path,
                self.requires()['taxonomy_file'].path, self.output_tree]
        with open(f'{self.output_tree}-log', 'w') as fh:
            proc = subprocess.Popen(args, stdout=fh, encoding='utf-8')
            proc.communicate()

    def output(self):
        return {'tree': luigi.LocalTarget(self.output_tree),
                'table': luigi.LocalTarget(f'{self.output_tree}-table'),
                'tax': luigi.LocalTarget(f'{self.output_tree}-taxonomy'),
                'log': luigi.LocalTarget(f'{self.output_tree}-log')}

    def get_tax(self):
        out = dict()
        with self.output()['tax'].open('r') as fh:
            for line in fh.readlines():
                gid, tax = line.strip().split('\t')
                out[gid] = tax.replace('; ', ';')
        return out


class GenerateCommonTaxa(luigi.Task):
    """Generates a list of common taxa for all input trees.
    GenerateCommonTaxa(dir_out=self.dir_out, path_trees=self.path_trees)
    """
    path_trees = luigi.Parameter()
    dir_out = luigi.Parameter()

    def requires(self):
        return InputTreesFile(self.path_trees)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out, 'intermediate_results/common_taxa.tsv'))

    def common_taxa(self):
        common = set()
        with self.output().open('r') as fh:
            for line in fh.readlines():
                common.add(line.strip())
        return common

    def run(self):
        common_taxa = None
        for name, path in self.requires().dict_name_path().items():
            tree = dendropy.Tree.get_from_path(path,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)
            taxa = {x.label for x in tree.taxon_namespace}
            if common_taxa:
                common_taxa = common_taxa.intersection(taxa)
            else:
                common_taxa = taxa

        with self.output().open('w') as fh:
            for taxon in sorted(common_taxa):
                fh.write(f'{taxon}\n')


class TreeDistances(luigi.Task):
    """Compute the robinson-foulds distance metric between two trees
        return TreeDistances(dir_out=self.dir_out,
                             tree_1_id='00.ar122_r89_reference',
                             tree_1_path='/home/uqamussi_local/github/treecompare/example/trees/00.ar122_r89_reference_rooted.tree',
                             tree_2_id='14.1.SSU.FastTree.900bp',
                             tree_2_path='/home/uqamussi_local/github/treecompare/example/trees/14.1.SSU.FastTree.900bp_rooted.tree',
                             common='global',
                             path_trees=self.path_trees)')
                             """
    tree_1_id = luigi.Parameter()
    tree_2_id = luigi.Parameter()
    tree_1_path = luigi.Parameter()
    tree_2_path = luigi.Parameter()
    common = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_trees = luigi.Parameter(default=None)  # must be specified with common

    def requires(self):
        out = {'tree_1': NewickTreeFile(path=self.tree_1_path),
               'tree_2': NewickTreeFile(path=self.tree_2_path)}
        if self.common == 'global':
            out['common_taxa'] = GenerateCommonTaxa(dir_out=self.dir_out, path_trees=self.path_trees)
        return out

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out,
                                              'intermediate_results',
                                              f'robinson_foulds_{self.common}',
                                              f'{self.tree_1_id}_vs_{self.tree_2_id}.tsv'))

    def run(self):

        genetreetk = GeneTreeTk()
        if self.common == 'global':
            rf, normalized_rf, num_taxa = genetreetk.robinson_foulds(self.requires()['tree_1'].output().path,
                                                           self.requires()['tree_2'].output().path,
                                                           self.requires()['common_taxa'].output().path)
        else:
            rf, normalized_rf, num_taxa = genetreetk.robinson_foulds(self.requires()['tree_1'].output().path,
                                                           self.requires()['tree_2'].output().path,
                                                           None)

        with self.output().open('w') as fh:
            fh.write(f'rf\t{rf}\n')
            fh.write(f'norm_rf\t{normalized_rf}\n')
            fh.write(f'n_common_taxa\t{num_taxa}\n')

    def get_metric(self):
        out = dict()
        with self.output().open('r') as fh:
            for line in fh.readlines():
                k, v = line.strip().split('\t')
                out[k] = v
        return out


class AggregateTreeDistances(luigi.Task):
    common = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_trees = luigi.Parameter()

    def requires(self):
        input_tree_file = InputTreesFile(path=self.path_trees)
        tree_dict = input_tree_file.dict_name_path()
        tree_ids = sorted(tree_dict.keys())
        for i in range(0, len(tree_ids)):
            tree_1_id = tree_ids[i]
            tree_1_path = tree_dict[tree_1_id]
            for j in range(0, i):
                tree_2_id = tree_ids[j]
                tree_2_path = tree_dict[tree_2_id]
                yield TreeDistances(dir_out=self.dir_out,
                                    tree_1_id=tree_1_id,
                                    tree_1_path=tree_1_path,
                                    tree_2_id=tree_2_id,
                                    tree_2_path=tree_2_path,
                                    common=self.common,
                                    path_trees=self.path_trees)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out,
                                              'results',
                                              f'robinson_foulds_{self.common}_mat.mat'))

    def run(self):
        self.output().makedirs()
        metrics = defaultdict(dict)
        names = set()
        for tc in self.requires():
            metrics[tc.tree_1_id][tc.tree_2_id] = tc.get_metric()['norm_rf']
            names.add(tc.tree_1_id)
            names.add(tc.tree_2_id)

        labels = sorted(list(names))
        mat = np.zeros([len(labels), len(labels)])
        for i in range(len(labels)):
            for j in range(i):
                if labels[i] in metrics and labels[j] in metrics[labels[i]]:
                    mat[i][j] = metrics[labels[i]][labels[j]]
        mat = mat + mat.T

        with self.output().open('w') as fh:
            fh.write('\t'.join(labels) + '\n')
            for i in range(len(labels)):
                fh.write('\t'.join(map(str, mat[i])) + '\n')

    def as_matrix(self):
        with self.output().open('r') as fh:
            labels = fh.readline().strip().split('\t')
            mat = np.zeros([len(labels), len(labels)])
            for i in range(len(labels)):
                row = fh.readline().strip().split('\t')
                for j in range(len(labels)):
                    mat[i][j] = float(row[j])
        return labels, mat


class CreateTreeFromTreeDistances(luigi.Task):
    """return CreateTreeFromTreeDistances(dir_out=self.dir_out, path_trees=self.path_trees, common='global')"""
    common = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_trees = luigi.Parameter()

    def requires(self):
        return AggregateTreeDistances(common=self.common, dir_out=self.dir_out,
                                      path_trees=self.path_trees)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out,
                                              'results',
                                              f'robinson_foulds_{self.common}.tree'))

    def run(self):
        self.output().makedirs()
        labels, mat = self.requires().as_matrix()

        # Convert the numpy distance matrix to the expected format
        mat_list = list()
        for i in range(len(labels)):
            row = list()
            for j in range(i + 1):
                row.append(mat[i][j])
            mat_list.append(row)

        dm = DistanceMatrix(names=labels, matrix=mat_list)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        # Write the tree to disk.
        with self.output().open('w') as fh:
            Phylo.write(tree, fh, 'newick')


class CreateClustermapFromTreeDistances(luigi.Task):
    """
    CreateClustermapFromTreeDistances(dir_out=self.dir_out, path_trees=self.path_trees, common='global')

    """
    common = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_trees = luigi.Parameter()

    def requires(self):
        return AggregateTreeDistances(common=self.common, dir_out=self.dir_out,
                                      path_trees=self.path_trees)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out,
                                              'results',
                                              f'robinson_foulds_{self.common}.svg'))

    def run(self):
        self.output().makedirs()
        labels, mat = self.requires().as_matrix()
        cmap = sns.cubehelix_palette(100, reverse=True)

        sns.set(font_scale=1)
        fig_size = (15, 15)

        rf_df = pd.DataFrame(mat, columns=labels, index=labels)
        sns.clustermap(rf_df, annot=True, fmt='.3f', square=True, cmap=cmap, figsize=fig_size).fig.suptitle(
            'Normalised Robinson–Foulds Distance')

        with self.output().open('w') as fh:
            plt.savefig(fh, format="svg")


class CompareNumOfTaxonDifferent(luigi.Task):
    """        return CompareNumOfTaxonDifferent(dir_out=self.dir_out, path_trees=self.path_trees,
                                          path_taxonomy=self.path_taxonomy)
"""
    path_trees = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_taxonomy = luigi.Parameter()

    def requires(self):
        out = {'true_tax': TaxonomyFile(self.path_taxonomy),
               'test_tax': list()}
        for name, path in InputTreesFile(path=self.path_trees).dict_name_path().items():
            out['test_tax'].append(DecorateTree(
                tree_id=name,
                input_tree=path,
                taxonomy_file=self.path_taxonomy,
                output_tree=os.path.join(self.dir_out, 'intermediate_results',
                                         'decorate', f'{name}.tree')))
        return out

    def output(self):
        return luigi.LocalTarget(os.path.join(self.dir_out, 'results',
                                              'number_of_taxon_different.tsv'))

    def run(self):
        self.output().makedirs()

        n_disagree, n_agree = dict(), defaultdict(lambda: 0)
        true_taxonomy = self.requires()['true_tax'].get_dict()

        for cur_model in self.requires()['test_tax']:
            n_disagree[cur_model.tree_id] = {'d': 0, 'p': 0, 'c': 0, 'o': 0,
                                             'f': 0, 'g': 0, 's': 0}
            test_taxonomy = cur_model.get_tax()

            for gid, test_tax in test_taxonomy.items():
                true_tax = true_taxonomy[gid]

                # Case 1: Agreement
                if true_tax == test_tax:
                    n_agree[cur_model.tree_id] += 1
                    continue

                # Case 2: Disagreement at a rank.
                true_ranks = convert_tax_to_ranks(true_tax)
                test_ranks = convert_tax_to_ranks(test_tax)
                for rank in ['d', 'p', 'c', 'o', 'f', 'g', 's']:

                    test_rank = test_ranks.get(rank)
                    true_rank = true_ranks.get(rank)

                    # Matching ranks.
                    if test_rank == true_rank:
                        continue

                    else:
                        n_disagree[cur_model.tree_id][rank] += 1
                        break

        # Write to disk
        with self.output().open('w') as fh:
            header = ['tree', 'domain', 'phylum', 'class', 'order', 'family', 'genus',
                      'species', 'total mismatch', 'total genomes', 'pct mismatch']
            fh.write('\t'.join(header) + '\n')
            for model_name, d_counts in n_disagree.items():
                tot_mm = sum([d_counts[x] for x in ('d', 'p', 'c', 'o', 'f', 'g', 's')])
                tot_g = n_agree[model_name] + tot_mm
                pct_mm = 100 * (tot_mm / tot_g)
                fh.write(f'{model_name}\t')
                fh.write('\t'.join([str(d_counts[x]) for x in ('d', 'p', 'c', 'o', 'f', 'g', 's')]) + '\t')
                fh.write(f'{tot_mm}\t{tot_g}\t{pct_mm:.2f}\n')


class CompareTreesPipeline(luigi.Task):
    path_trees = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_taxonomy = luigi.Parameter()
    common = luigi.Parameter()

    def requires(self):
        return CreateClustermapFromTreeDistances(dir_out=self.dir_out, path_trees=self.path_trees, common=self.common), \
               CompareNumOfTaxonDifferent(dir_out=self.dir_out, path_trees=self.path_trees,
                                          path_taxonomy=self.path_taxonomy), \
               CreateTreeFromTreeDistances(dir_out=self.dir_out, path_trees=self.path_trees, common=self.common)
