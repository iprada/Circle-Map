from setuptools import setup,find_packages

setup(name='Circle-Map',
      version='1.0',
      description='Circular DNA analysis tools',
      author='Inigo Prada-Luengo',
      url='https://github.com/iprada/Circle-Map',
      packages=find_packages(),
      headers=['/usr/include/zlib.h'],
      install_requires=[
          'pysam','pandas','biopython','numpy','pybedtools','edlib','numba','tqdm'
      ],
      entry_points={
          'console_scripts': [
              'Circle-Map = src.circle_map:main'
          ],

      },
     )
