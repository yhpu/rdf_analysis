from setuptools import setup, find_packages

setup(
    name="rdf_analysis",
    version="1.0.0",
    description="Radial Distribution Function Analysis Tool",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "matplotlib>=3.5.0",
        "scipy>=1.7.0",
        "MDAnalysis>=2.0.0"
    ],
    entry_points={
        'console_scripts': [
            'rdf-analysis=rdf_analysis.main:main'
        ]
    },
    python_requires=">=3.8",
)