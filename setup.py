from setuptools import setup, find_packages

setup(
    name='eatpy',
    version='0.1',
    author='Bolding-Bruggeman ApS',
    author_email='karsten@bolding-bruggeman.com',
    license='GPL',
    packages=find_packages(include=['eatpy*']),
    package_data={'eatpy': ['*.so', '*.dll', '*.dylib', '*.pyd']},
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'eat-filter = eatpy.filter:main',
            'eat-gotm-obs = eatpy.gotm.obs:main',
            'eat-gotm-gen = eatpy.gotm.generate_ensemble:main',
        ],
    }
)


