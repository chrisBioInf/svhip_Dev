
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="svhip_package_dev", 
    version="0.0.1",
    author="Christopher",
    author_email="christopher@bioinf.uni-leipzig.de",
    description="SVM training software for RNAz in development.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="TO UPDATE",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
