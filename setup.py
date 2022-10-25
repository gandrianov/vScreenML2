#!/usr/bin/env python
from setuptools import setup, find_packages

__version__ = "2.0.0"

setup(
    name="vScreenML",
    version=__version__,
    description="ML classifier for rescoring of virtual screening hits to prune out false positives",
    author="Grigorii Andrianov, Yusuf Adeshina and John Karanicolas",
    author_email="grigorii.andrianov@gmail.com",
    url="https://github.com/gandrianov/vScreenML",
    license="MIT",
    packages=find_packages(),
    package_data={"":["utils/*.txt",
                      "models/*",
                      "data/*",
                      "external/binana/*.md"]},
    include_package_data=True,
    install_requires=open('requirements.txt', 'r').readlines(),
    entry_points="""
        [console_scripts]
        vscreenml_calculate_features=vScreenML.main:calculate_features
        vscreenml_minimize_complex=vScreenML.main:minimize_pdb_complex
        vscreenml_mol2params=vScreenML.utils.mol2genparams_utils:run_mol2params
        vscreenml_predict_score=vScreenML.utils.xgboost_utils:predict_vscreenml_score
        vscreenml_train_model=vScreenML.utils.xgboost_utils:train_model
        """,
    keywords=['cheminformatics', 'virtual screening'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
