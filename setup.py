
from setuptools import setup

setup(name='chemkin',
      version='1.3',
      description='The best chemkin package',
      url='https://github.com/G12-cs207-FinalProject/cs207-FinalProject.git',
      author='Parser',
      license='Harvard',
      packages=['chemkin',
                'chemkin.preprocessing',
                'chemkin.preprocessing.tests',
                'chemkin.reaction',
                'chemkin.reaction.tests',
                'chemkin.thermodynamics'],
    package_data={'chemkin': ['*.sqlite',
                                'xml-files/*.xml']})