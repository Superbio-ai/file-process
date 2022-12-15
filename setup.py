from setuptools import setup

setup(name="file_process",
      version="1.0.3",
      description="A package that does file validation and file preview.",
      license='MIT',
      author="superbio.ai",
      url='https://github.com/Superbio-ai/file-process',
      install_requires=['pandas==1.2.5', 'anndata==0.8.0', 'Werkzeug==2.0.3', 'scanpy==1.9.1', 'numpy==1.21'],
      packages=['file_process'])
