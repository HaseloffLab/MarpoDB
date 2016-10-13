from setuptools import setup, find_packages
from codecs import open
from os import path


setup(
	name = 'PartsDB',
	version = '0.0.1.dev1',
	description='Database of genetic parts',
	url='https://github.com/HaseloffLab/MarpoDB',
	author = 'Mihails Delmans',
	author_email='md656@cam.ac.uk',
	license = 'GPL',
	classifiers=[
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 3 - Alpha',

		# Indicate who your project is intended for
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',

		# Pick your license as you wish (should match "license" above)
		 'License :: OSI Approved :: GNU General Public License (GPL)',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 2.7',
	],
	keywords='bioinfomratics parts synthetic biology',
	packages=find_packages(),
	install_requires=['sqlalchemy'],
)

