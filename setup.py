import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mat_tools",
    version="0.1",
    author="F. Millour, P. Berio, A. Meilland, J. Varga",
    author_email="florentin.millour@oca.eu",
    description="A package to help reducing MATISSE/VLTI data in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.oca.eu/MATISSE/tools",
    packages=['mat_tools'],
    package_data={'mat_tools':['icons/*']},
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=['objectlistview','astropy','wxpython','tqdm']
)
