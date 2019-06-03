#MIT License
#
#Copyright (c) 2019 Iñigo Prada Luengo
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

from setuptools import setup,find_packages




setup(name='Circle-Map',
      version='1.0.2',
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
              'Circle-Map = circlemap.circle_map:main'
          ],

      },
      classifiers=[
          'License :: OSI Approved :: MIT License'
      ],
     )
