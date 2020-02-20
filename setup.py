from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="meta_analysis",
    version="0.0.1",
    author="Kang Xiaopeng",
    description="Voxelwise Meta analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages('./meta_analysis'),
    python_requires='>=3.7',
)