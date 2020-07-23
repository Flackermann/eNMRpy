from setuptools import setup

setup(
   name='eNMRpy',
   version='0.0.2',
   author='Florian Ackermann',
   author_email='nairolf.ackermann@gmail.com',
   packages=['eNMRpy', 'eNMRpy.Measurement'],
   #scripts=['bin/script1','bin/script2'],
   #url='http://pypi.python.org/pypi/PackageName/',
   license='LICENSE.txt',
   description='nmrglue-based package for the import and analysis of electrophoretic NMR-data',
   long_description=open('README.txt').read(),
   install_requires=[
       "nmrglue",
       "lmfit",
       "pandas",
       "matplotlib",
       "numpy",
       "scipy",
       'scikit-learn<=0.23',
   ],
)
