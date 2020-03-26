from setuptools import setup

setup(name='treecompare',
      version='0.1.0',
      packages=['treecompare'],
      entry_points={
          'console_scripts': [
              'treecompare = treecompare.__main__:main'
          ]
      },
      install_requires=['pandas', 'luigi', 'dendropy', 'ete3', 'matplotlib', 'numpy', 'biopython', 'seaborn'],
      )
