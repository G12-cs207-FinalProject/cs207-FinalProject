from setuptools import setup

AUTHOR = 'P&P'
VERSION = 0.5

setup(name='chemkin',
      version=VERSION,
      description='The best chemkin package',
      url='https://github.com/G12-cs207-FinalProject/cs207-FinalProject.git',
      author=AUTHOR,
      license='Harvard',
      packages=['chemkin',
                'chemkin.preprocessing',
                'chemkin.preprocessing.tests',
                'chemkin.reaction',
                'chemkin.reaction.tests',
                'chemkin.thermodynamics'],
      package_data={'chemkin': ['*.sqlite',
                                'xml-files/*.xml']
                    }
      )
