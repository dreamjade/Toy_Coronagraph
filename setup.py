from setuptools import setup, find_packages, Extension
import re

# auto-updating version code stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

def get_requires():
    reqs = []
    for line in open("requirements.txt", "r").readlines():
        reqs.append(line)
    return reqs

setup(
    name="toycoronagraph",
    version=get_property("__version__", "toycoronagraph"),
    description="Toy!!",
    url="https://github.com/dreamjade/Toy_Coronagraph",
    author="dreamjade",
    author_email="",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
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
        "Programming Language :: Python :: 3.11",
    ],
    install_requires=get_requires(),
    keywords="Toy Coronagraph"
)
