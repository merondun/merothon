from setuptools import setup, find_packages

setup(
    name='merothon',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'calculate_r2=merothon.Calculate_R2_LD_VCFs:main',
            'plot_genotypes=merothon.Plot_Genotypes:main',
            'permutation_test=merothon.Permutation_Test:main',
            'bootstrap_sample=merothon.Bootstrap_Sample_Region:main',
        ],
    },
    install_requires=[
        'pandas',
        'pysam',
        'matplotlib',
	'scipy',
	'seaborn'
        # Add other dependencies as needed
    ],
)

