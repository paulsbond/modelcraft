import setuptools
from setuptools.command.install import install
import os

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()


class InstallPackage(install):
    def __init__(self, dist):
        super(install, self).__init__(dist)

    def run(self):
        install.run(self)
        self.__post_install()

    def __post_install(self):
        os.system('nucleofind-install --all')


setuptools.setup(
    name="modelcraft",
    version="5.0.0",
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
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires="~=3.7",
    install_requires=[
        "gemmi >=0.5.4",
        "numpy",
        "pandas",
        "requests",
        "scipy",
        "nucleofind"
    ],
    entry_points={
        "console_scripts": [
            "modelcraft = modelcraft.scripts.modelcraft:main",
            "modelcraft-contents = modelcraft.scripts.contents:main",
            "modelcraft-copies = modelcraft.scripts.copies:main",
        ]
    },
    cmdclass={'install': InstallPackage}
)
