from setuptools import setup

setup(name='chemkin',
      version='3.4',
      description='The best chemkin package',
      url='https://github.com/G12-cs207-FinalProject/cs207-FinalProject.git',
      author='Monsieur Package du Parser',
      author_email='nathaniel_stein@g.harvard.edu',
      license='Harvard',
      packages=['chemkin',
                'chemkin.preprocessing',
                'chemkin.preprocessing.tests',
                'chemkin.reaction',
                'chemkin.reaction.tests',
                'chemkin.thermodynamics',
                'chemkin.thermodynamics.tests',
                'chemkin.viz',
                'chemkin.viz.tests',
                'chemkin.solver',
                'chemkin.solver.tests'],
      package_data={'chemkin':['thermodynamics/*.sqlite',
                               'xml-files/*.xml']})

