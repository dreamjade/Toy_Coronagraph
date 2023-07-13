from setuptools import setup, find_packages
setup(
    name="toycoronagraph",
    version=get_property("__version__", "toycoronagraph"),
    description="Toy!!",
    url="https://github.com/dreamjade/Toy_Coronagraph",
    author="dreamjade",
    author_email="",
    license="MIT",
    packages=find_packages(),
    #package_data={"": ["kernels/*.cu"]},
    #ext_modules=get_extensions(),
    #include_dirs=[numpy.get_include()],
    #include_package_data=True,
    zip_safe=False,
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
    ],
    #install_requires=get_requires(),
    keywords="Toy Coronagraph"
)