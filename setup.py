from setuptools import setup

setup(name='truss',
      version='0.1',
      description="Scriptable 2D truss framework for teaching purposes",
      long_description="",
      author='Ludovic Charleux',
      author_email='ludovic.charleux@univ-smb.fr',
      license='GPL v2',
      packages=['truss'],
      zip_safe=False,
      install_requires=[
          "numpy",
          "scipy",
          "matplotlib"
          ],
      )
