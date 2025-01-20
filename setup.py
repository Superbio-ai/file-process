from setuptools import setup
import pathlib

HERE = pathlib.Path(__file__).parent

README = (HERE / "README.md").read_text()

setup(name="file_process",
      version="1.4.5",
      description="A package that does file validation and file preview.",
      long_description=README,
      long_description_content_type="text/markdown",
      license='MIT',
      author="superbio.ai",
      url='https://github.com/Superbio-ai/file-process',
      install_requires=['pandas==2.2', 'anndata==0.10.5', 'numpy>=1.21'],
      packages=['file_process', 'file_process.h5ad', 'file_process.csv', 'file_process.txt'])
