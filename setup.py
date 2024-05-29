from setuptools import setup, find_packages

setup(
    name="MicEASH",
    version="0.1.0",
    author="Hiroaki Ohishi",
    author_email="hirohishi@outlook.jp",
    description="A command line tool for DNA-seqFISH+ analysis",
    url="https://github.com/hirohishi/MicEASH",  # Optional: Your project's main homepage.
    packages=find_packages(),  # Automatically find and include all packages
    install_requires=[
        "pandas>=1.0.0",
        "numpy>=1.18.0",
        "scipy>=1.4.0",
        "seaborn>=0.10.0",
        "matplotlib>=3.1.0",
        "scikit-learn>=0.22.0",
        "click>=7.0",
        "openpyxl>=3.0.0",
        "ipykernel>=5.1.0",

  # Ensure you include all necessary dependencies
    ],
    entry_points={
        "console_scripts": [
            "MicEASH=mic_eash.cli:cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

