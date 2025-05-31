from setuptools import setup, find_packages

setup(
    name="isoqmap",
    version="0.1.0",
    description="Your project description",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/IsoQMap",  # 如果有github
    license="MIT",  # 或其他许可证
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=[
        # 填写依赖包，比如 "click", "numpy" 等
    ],
    entry_points={
        "console_scripts": [
            "isoqmap=isoqmap.main:main",  # 命令行入口
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
