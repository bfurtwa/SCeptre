from setuptools import setup, find_packages

setup(
    name='sceptre',
    version='0.1',
    description='Analysis of Multiplexed Single-Cell Proteomics Data',
    url='https://github.com/bfurtwa/sceptre',
    download_url='https://github.com/bfurtwa/sceptre/archive/v0.1.tar.gz',
    packages=find_packages(exclude=['Schoof_et_al', '.idea']),
    python_requires='>=3.6',
    install_requires=[
        'scanpy',
        'typing_extensions; python_version < "3.8"',  # for `Literal`,
        'adjustText',
    ],
    author='Benjamin Furtwängler',
    author_email='bfurtwaengler@web.de',
    license='MIT'
)