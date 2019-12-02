import setuptools

with open("README.md", "r") as f:
  long_description = f.read()

setuptools.setup(
  name="autocoord",
  version="0.0.0",
  author="Paul Bond",
  author_email="paul.bond@york.ac.uk",
  description="Automated model building pipeline",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/paulsbond/autocoord",
  packages=setuptools.find_packages(),
  classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)"
    "Operating System :: OS Independent",
  ],
  python_requires=">=3.6",
)
