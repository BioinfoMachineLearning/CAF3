import copy
import ml_collections
MONOMER_PREDICTIONS_PER_MODEL = 100



af3_program_path = "/bmlfast/pngkg/Benchmark_AF3/alphafold3"
af3_params_path = "/home/pngkg/public_databases/params/"
af3_db_path = "/home/pngkg/public_databases/"
env_dir = "/bmlfast/bml_casp16/anaconda3/envs/multicom4/bin/"
database_dir = "/bmlfast/bml_casp16/tools/alphafold_databases_multicom3"

pdb_mmcif_database_dir = f"{database_dir}/pdb_mmcif//mmcif_files/"

MONOMER_CONFIG = ml_collections.ConfigDict({
    'common_config': {
        'msa_source': 'default',
        'template_source': 'pdb70',
    },
    'predictors':{
        'default': {
        },
        'default_pdb70_new': {
            'template_source': 'pdb70_newest',
        },
        'default_seq_temp': {
            'template_source': 'pdb_sort90'
        },
        'original': {
            'msa_source': 'original',
        },
        'ori_seq_temp': {
            'msa_source': 'original',
            'template_source': 'pdb_sort90'
        },
        'colabfold_web': {
            'msa_source': 'colabfold_web',
        },
        'colabfold_web_not': {
            'msa_source': 'colabfold_web',
            'template_source': 'notemplate',
        },
        'dhr': {
            'msa_source': 'dhr',
        },
        'def_notemp': {
            'template_source': 'notemplate',
        },
        'deepmsa_dMSA_hhb':{
            'msa_source': 'dMSA.hhb',
        },
        'deepmsa_dMSA_jac':{
            'msa_source': 'dMSA.jac',
        },
        'deepmsa_dMSA_hms':{
            'msa_source': 'dMSA.hms',
        },
        'deepmsa_dMSA':{
            'msa_source': 'dMSA',
        },
        'deepmsa_qMSA':{
            'msa_source': 'qMSA',
        },
        'deepmsa_aMSA':{
            'msa_source': 'aMSA',
        },
        'deepmsa_qMSA_hhb':{
            'msa_source': 'qMSA.hhb',
        },
        'deepmsa_qMSA_jac':{
            'msa_source': 'qMSA.jac',
        },
        'deepmsa_qMSA_hh3':{
            'msa_source': 'qMSA.hh3',
        },
        'deepmsa_qMSA_hms':{
            'msa_source': 'qMSA.hms',
        },
        'deepmsa_DeepJGI_hms':{
            'msa_source': 'DeepJGI.hms',
        },
        'deepmsa_DeepJGI':{
            'msa_source': 'DeepJGI',
        },
        'deepmsa_q3JGI':{
            'msa_source': 'q3JGI',
        },
        'deepmsa_q4JGI':{
            'msa_source': 'q4JGI',
        },
        'deepmsa_q3JGI_hms':{
            'msa_source': 'q3JGI.hms',
        },
        'deepmsa_q4JGI_hms':{
            'msa_source': 'q4JGI.hms',
        },
        'def_esm_msa': {
            'input_msa_source': 'default',
            'msa_source': 'esm_msa',
            'num_ensemble': 1,
            'num_recycle': 3,
        },
        'def_dom_hhsearch':{
            'start_msa': 'default',
            'domain_msa_source': 'default',
            'msa_source': 'def_dom_hhsearch',
        },
        'def_dom_parser':{
            'start_msa': 'default',
            'domain_msa_source': 'default',
            'msa_source': 'def_dom_parser',
        },
        'def_dom_unidoc':{
            'start_msa': 'default',
            'domain_msa_source': 'default',
            'msa_source': 'def_dom_unidoc',
        },
        'def_dom_manual':{
            'start_msa': 'default',
            'domain_msa_source': 'default',
            'msa_source': 'def_dom_manual',
        },
        'dmsa_dom_hhsearch':{
            'start_msa': 'deepmsa2_dMSA',
            'domain_msa_source': 'deepmsa2_dMSA',
            'msa_source': 'dmsa_dom_hhsearch',
        },
        'dmsa_dom_parser':{
            'start_msa': 'deepmsa2_dMSA',
            'domain_msa_source': 'deepmsa2_dMSA',
            'msa_source': 'dmsa_dom_parser',
        },
        'dmsa_dom_unidoc':{
            'start_msa': 'deepmsa2_dMSA',
            'domain_msa_source': 'deepmsa2_dMSA',
            'msa_source': 'dmsa_dom_unidoc',
        },
        'dmsa_dom_manual':{
            'start_msa': 'deepmsa2_dMSA',
            'domain_msa_source': 'deepmsa2_dMSA',
            'msa_source': 'dmsa_dom_manual',
        },
    }
})
