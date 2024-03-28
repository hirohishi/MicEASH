from setuptools import setup

setup(
    name="MicEASH",
    version="0.0.1",
    author="Hiroaki Ohishi",
    author_email = "hirohishi@outlook.jp",
    install_requires=[],
    entry_points={
        "console_scripts": [
            "MicEASH = src:cli"
        ]
    }
)
