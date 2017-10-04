from setuptools import setup

setup(
    name='tidyss',
    description='Tools for better sample sheets',
    version='0.1',
    author='Ryan Richholt',
    author_email='rrichholt@tgen.org',
    url='https://github.com/ryanrichholt/tidyss',
    packages=['tidyss'],
    entry_points={
        'console_scripts': [
            'tidyss = tidyss.fastqs:main'
        ]
    },
    install_requires=[
        'pyyaml',
    ],
)