from setuptools import setup

setup(
   name='eNMRpy',
   version='0.0.6',
   author='Florian Ackermann',
   author_email='nairolf.ackermann@gmail.com',
   packages=['eNMRpy', 'eNMRpy.Measurement'],
   url='https://github.com/Flackermann/eNMRpy/',
   license='LICENSE.txt',
   description='nmrglue-based package for the import and analysis of electrophoretic NMR-data',
   long_description=open('README.md').read(),
   install_requires=[
       "nmrglue",
       "lmfit>=1.0.1",
       "pandas",
       "matplotlib>=3.1.3",
       "numpy",
       "scipy",
       'scikit-learn<=0.23',
   ],
)
