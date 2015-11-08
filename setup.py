from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='TwoDimBoussinesq',
      version='0.1dev1',
      description='A pythonian solver for 2D Boussinesq equations\
              analysis',
      url='.',
      author='The pyboussinesq team',
      author_email='',
      license='MIT',
      packages=['TwoDimBoussinesq'],
      install_requires=[
          'numpy',
      ],
      test_suite = 'nose.collector',
      zip_safe=False)
