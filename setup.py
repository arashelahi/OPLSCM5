from setuptools import setup, find_packages

setup(
    name="OPLSCM5",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",   # often installed as 'openbabel-wheel' or from conda
    ],
    entry_points={
        "console_scripts": [
            "itp_gen = OPLSCM5.itp_rewrite_main:main"
        ]
    }
)
