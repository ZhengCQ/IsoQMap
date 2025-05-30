from setuptools import setup, find_packages
import os

# 自动检测是否为Conda环境
is_conda = 'CONDA_PREFIX' in os.environ

data_files = []
if is_conda:
    # Conda环境安装到envs内的share目录
    conda_share = os.path.join(os.environ['CONDA_PREFIX'], 'share', 'isoqmap')
    data_files.append((conda_share, ['data/config/default.ini']))
else:
    # 常规Linux系统安装
    data_files.append(('/usr/local/share/isoqmap', ['data/config/default.ini']))


setup(
    name='isoqmap',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='Splicing-regulatory Driver Genes Identification Tool',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/asge-demo',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'click',
        'pandas',
        'numpy',
        'scipy',
        'statsmodels',
    ],
    entry_points={
        'console_scripts': [
            'isoqmap=isoqmap.main:cli'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
