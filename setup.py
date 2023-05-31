from setuptools import setup, find_packages

setup(
    name='eatpy',
    version='0.9.4',
    author='Bolding-Bruggeman ApS',
    author_email='karsten@bolding-bruggeman.com',
    license='GPL',
    packages=find_packages(include=['eatpy*']),
    package_data={'eatpy': ['*.so', '*.dll', '*.dylib', '*.pyd']},
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'eat-gotm-gen = eatpy.models.gotm_generate_ensemble:main',
        ],
    }
)


