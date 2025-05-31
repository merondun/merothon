from setuptools import setup, find_packages

setup(
    name='half_merothon',
    version='0.4.2',
    packages=find_packages(),
    description='Genomic utility scripts',
    author='J Merondun',
    entry_points={
        'console_scripts': [
            'vcf_to_pca=merothon.VCF_to_PCA:main',
            'plot_ld=merothon.Plot_LD:main',
            'calculate_r2=merothon.Calculate_R2_LD_VCFs:main',
            'plot_genotypes=merothon.Plot_Genotypes:main',
            'permutation_test=merothon.Permutation_Test:main',
            'polarize_vcf=merothon.Assign_Ancestral_Allele:main',
            'count_mutations=merothon.Calculate_Fasta_Mutations:main',
            'map_chromosomes=merothon.Reference_to_Scaffold_ChrID:main',
            'subset_snps=merothon.Subset_SNPs:main',
        ],
    },
    install_requires=[
        'pandas',
        'numpy',
        'pysam',
        'matplotlib-base',
        'biopython',
#        'scikit-allel',
    ],
)

