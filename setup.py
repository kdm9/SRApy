from setuptools import setup
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = 'srapy/_version.py'
versioneer.versionfile_build = 'srapy/_version.py'
versioneer.tag_prefix = ''
versioneer.parentdir_prefix = 'srapy-'

desc = """
SRApy: Scripts to download SRA files
"""

install_requires = [
    "biopython==1.65",
    "docopt==0.6.2",
    "lxml==3.4.4",
    "progressbar2==2.7.3",
]

test_requires = [
    "coverage==3.7.1",
    "nose==1.3.0",
    "pep8==1.4.6",
    "pylint==1.0.0",
]

setup(
    name="srapy",
    packages=['srapy', ],
    scripts=['scripts/get-project-sras.py',],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    install_requires=install_requires,
    tests_require=test_requires,
    description=desc,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/srapy",
    keywords=[
        "bioinformatics",
        "SRA",
        "Next-gen Sequencing",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later " +
            "(GPLv3+)",
    ],
    test_suite="test",
)
