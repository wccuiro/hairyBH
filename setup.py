from setuptools import setup

readme = open("README.md", "r")


setup(
    name="hairyBH",
    packages=['hairyBH'],
    version="0.1.1",
    description="A package to calculate some important parameters of hairy black holes",
    long_description=readme.read(),
    long_description_content_type="text/markdown",
    author="Weyner Ccuiro",
    author_email="wey.090300@gmail.com",
    url="https://github.com/Skylines316/hairyBH",
    download_url="https://github.com/Skylines316/hairyBH/tarball/0.1.1",
    keywords=['blackholes','hairy','geodesics'],
    classifiers=[],
    license="MIT",
    include_package_data=True,
    install_requires=[
        "numpy>=1.21",
        "scipy>=1.7",
    ],
)