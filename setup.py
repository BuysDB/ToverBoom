from setuptools import setup
setup(
	name='toverboom',
	version='0.01',
	author='Buys de Barbanson',
	author_email='code@buysdb.nl',
	description='A python module to generate good looking lineage trees ',
	url='https://github.com/BuysDB/ToverBoom',
	py_modules=['toverboom'],
	install_requires=['matplotlib','numpy','networkx','pandas']
)
