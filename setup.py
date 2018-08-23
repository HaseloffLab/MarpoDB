from setuptools import setup, find_packages
from codecs import open
from os import path

setup(
	name = 'marpodb',
	version = '4.0',
	description='Marchantia polymorpha gene-centric database',
	# url='https://github.com/HaseloffLab/marpodb',
	# download_url='https://github.com/HaseloffLab/LoopDB/archive/0.2.3.tar.gz',
	author = 'Mihails Delmans',
	author_email='md656@cam.ac.uk',
	license = 'MIT',
	classifiers=[
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 4 - Beta',

		# Indicate who your project is intended for
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',

		# Pick your license as you wish (should match "license" above)
		 'License :: OSI Approved :: MIT License',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 2.7',
	],
	keywords='bioinfomratics',
	packages=find_packages(),
	install_requires=[
		'sqlalchemy', 'partsdb',
		'biopython', 'psycopg2-binary',
		'wtforms', 'flask', 'flask-user',
		'Flask-SQLAlchemy', 'flask_bcrypt', 'Flask-Restless', 'sqlalchemy-utils',
		'sh'
	],
	package_data={
			"marpodb":
			[
				"server/templates/*",
				"server/static/stylesheets/*",
				"server/static/javascript/*",
				"server/static/img/*",
				"server/static/json/*",
				"server/static/interpro/*",
				"server/static/interpro/resources/*",
				"server/static/interpro/resources/css/*",
				"server/static/interpro/resources/freeMarker/*",
				"server/static/interpro/resources/images/*",
				"server/static/interpro/resources/javascript/*",
				"server/static/interpro/resources/javascript/qtip2/*"
			]
	},
	include_package_data=True,
	scripts=['scripts/installmarpodb']
)

