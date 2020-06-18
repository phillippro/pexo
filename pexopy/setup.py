from setuptools import setup

setup(
    name="pexopy",
    url="https://github.com/timberhill/pexopy",
    author="Maksym Lisogorskyi",
    author_email="m.lisogorskyi@gmail.com",
    packages=["pexopy"],
    install_requires=["numpy"],
    python_requires='>=3.6',
    version="0.1",
    license="MIT",
    description="A python wrapper for PEXO software",
    long_description=open("README.md").read(),
    classifiers=[
        "Development Status :: 3 - Alpha"
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords="Barycentric Correction Astronomy Radial Velocity Astrometry Timing General Relativity",
)