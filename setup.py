from setuptools import setup,find_packages
import os



setup(name='Circle-Map',
      version='1.0',
      description='Circular DNA analysis tools',
      author='Inigo Prada-Luengo',
      url='https://github.com/iprada/Circle-Map',
      packages=find_packages(),
      install_requires=[
          'pysam>=0.15.2','pybedtools>=0.8.0','pandas>=0.24.2','biopython>=1.73','numpy>=1.16.3',
          'edlib>=1.2.4','numba>=0.43.1','tqdm>=4.31.1','scipy>=1.2.1'
      ],
      entry_points={
          'console_scripts': [
              'Circle-Map = src.circle_map:main'
          ],

      },
     )
