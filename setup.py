import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="prodes",
    version="1.0.0",
    author="Tim Neijenhuis",
    author_email="t.neijenhuis@tudelft.nl",
    description="A package to calculate protein surface descriptors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tneijenhuis/github",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    package_requires=[
        "numpy", "pandas"
    ],
    package_data={"prodes": [
        "data/*"
    ]}
)