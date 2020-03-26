import luigi


class InputTreesFile(luigi.ExternalTask):
    """A user-specified batch file specifying the tree id and trees."""
    path = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.path)

    def dict_name_path(self):
        trees = dict()
        with self.output().open('r') as fh:
            for line in fh.readlines():
                name, path = line.strip().split('\t')
                if name.startswith('#'):
                    continue
                trees[name] = path
        return trees


class NewickTreeFile(luigi.ExternalTask):
    """A user-specified newick representing at tree."""
    path = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.path)


class TaxonomyFile(luigi.ExternalTask):
    """A taxonomy file provided by the user as input."""
    path = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.path)

    def get_dict(self):
        out = dict()
        with self.output().open('r') as fh:
            for line in fh.readlines():
                gid, tax = line.strip().split('\t')
                out[gid] = tax.replace('; ', ';')
        return out
