import setuptools

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="modelcraft",
    version="0.2.8",
    author="Paul Bond",
    author_email="paul.bond@york.ac.uk",
    description="Automated model building pipeline for X-ray crystallography",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/paulsbond/modelcraft",
    packages=setuptools.find_packages(),
    include_package_data=True,
    license="LGPL-2.1",
    classifiers=[
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires=">=3.6",
    install_requires=["gemmi >=0.4.0", "pandas"],
    entry_points={"console_scripts": ["modelcraft = modelcraft.__main__:main"]},
)
