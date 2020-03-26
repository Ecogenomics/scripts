import os

import luigi

from treecompare.f_measure_tree import FMeasureTree
from treecompare.io import InputTreesFile
from treecompare.pipeline.compare_trees import DecorateTree


class GenerateFMeasureTree(luigi.Task):
    path_trees = luigi.Parameter()
    dir_out = luigi.Parameter()
    path_taxonomy = luigi.Parameter()
    plot_legend = luigi.BoolParameter(default=True)

    def requires(self):
        for name, path in InputTreesFile(path=self.path_trees).dict_name_path().items():
            yield DecorateTree(
                tree_id=name,
                input_tree=path,
                taxonomy_file=self.path_taxonomy,
                output_tree=os.path.join(self.dir_out, 'intermediate_results',
                                         'decorate', f'{name}.tree'))

    def output(self):
        f_hash = abs(hash(frozenset({x.tree_id for x in self.requires()})))
        return luigi.LocalTarget(os.path.join(self.dir_out, 'results',
                                              f'fmeasure_tree_{f_hash}.svg'))

    def run(self):
        self.output().makedirs()
        fmt = FMeasureTree(self.path_taxonomy)

        for model in self.requires():
            fmt.add_table(model.tree_id, model.output()['table'].path)

        fmt.run(legend=self.plot_legend, out_path=self.output().path)

