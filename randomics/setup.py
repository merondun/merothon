from setuptools import setup, find_packages

setup(
    name='merothon',
    version='0.4.2',
    packages=find_packages(),
    description='Genomic utility scripts',
    author='J Merondun',
    entry_points={
        'console_scripts': [
            'vcf_to_pca=merothon.scripts.VCF_to_PCA:main',
            'plot_ld=merothon.scripts.Plot_LD:main',
            'calculate_r2=merothon.scripts.Calculate_R2_LD_VCFs:main',
            'plot_genotypes=merothon.scripts.Plot_Genotypes:main',
            'permutation_test=merothon.scripts.Permutation_Test:main',
            'polarize_vcf=merothon.scripts.Assign_Ancestral_Allele:main',
            'count_mutations=merothon.scripts.Calculate_Fasta_Mutations:main',
            'map_chromosomes=merothon.scripts.Reference_to_Scaffold_ChrID:main',
            'subset_snps=merothon.scripts.Subset_SNPs:main',
        ],
    },
    install_requires=[
        'pandas',
        'numpy',
        'pysam',
        'matplotlib-base',
        'biopython',
        'scikit-allel',
    ],
)

