from setuptools import setup, find_packages
import sys
import codecs


def main():

    with open("requirements.txt") as f:
        requirements = f.read().splitlines()

    setup(
        maintainer="Rasmus Zetter",
        maintainer_email="rasmus.zetter@aalto.fi",
        name="bfieldtools",
        description="Open-source Python software for magnetic field modeling",
        long_description=codecs.open("README.rst", encoding="utf8").read(),
        url="https://bfieldtools.github.io/",
        license="BSD (3-clause)",
        download_url="https://github.com/bfieldtools/bfieldtools",
        use_scm_version={
            "write_to": "bfieldtools/version.py",
            "write_to_template": '__version__ = "{version}"',
        },
        setup_requires=["setuptools_scm"],
        python_requires=">=3.6",
        install_requires=requirements,
        packages=find_packages(),
        include_package_data=True,
        package_data={"bfieldtools": ["example_meshes/*"]},
        classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "Programming Language :: Python",
            "Topic :: Software Development",
            "Topic :: Scientific/Engineering",
        ],
        platforms="any",
        zip_safe=False,
    )


if __name__ == "__main__":
    main()
