from setuptools import find_packages, setup


setup(
    name="nddia",
    version=0.0.1,
    description="N/D diabatization in python.",
    url="https://github.com/eljost/nddia",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "matplotlib",
        "numpy",
        "pytest",
    ],
)
