import os
import sys
import json
import argparse
import gzip
import configparser
import pandas as pd

def main():

    Studies_to_run=['GCST006979', 'GCST90025948', 'GCST90038635', 'GCST90038637',
        'GCST004131', 'GCST90002298', 'GCST007430', 'GCST90019464',
        'GCST90012111', 'GCST90012110', 'GCST90012109', 'GCST90012108',
        'GCST90025958', 'GCST90002316', 'GCST90000026', 'GCST90000025',
        'GCST90000027', 'GCST90038595', 'GCST90002405', 'GCST90002406',
        'GCST90002310', 'GCST004607', 'GCST90002391', 'GCST90002328',
        'GCST90002397', 'GCST90002402', 'GCST004603', 'GCST90002357',
        'GCST90002369', 'GCST007429', 'GCST90013405', 'GCST90002382',
        'GCST90025988', 'GCST90002412', 'GCST90002400', 'GCST90002404',
        'GCST90002395', 'GCST90002346', 'GCST90002385', 'GCST90002386',
        'GCST007432', 'GCST90002388', 'GCST90002396', 'GCST90002363',
        'GCST004616', 'GCST90002401', 'GCST005194', 'GCST90002392',
        'GCST90002334', 'GCST90002390', 'GCST004630', 'GCST90002322',
        'GCST010144', 'GCST90000616', 'GCST90000614', 'GCST90000617',
        'GCST90000615', 'GCST90000618', 'GCST006941', 'GCST90038616',
        'GCST90002381', 'GCST005531', 'GCST005523', 'GCST006288',
        'GCST90002340', 'GCST90002394', 'GCST90002387', 'GCST004599',
        'GCST004602', 'GCST90002389', 'GCST90002403', 'GCST004601',
        'GCST90038638', 'GCST90038646', 'GCST90002304', 'GCST90002384',
        'GCST90002383', 'GCST004604', 'GCST004615', 'GCST90000059',
        'GCST007431', 'GCST004612', 'GCST004619', 'GCST004621',
        'GCST006804', 'GCST004628', 'GCST004605', 'GCST004611',
        'GCST004622', 'GCST006950', 'GCST90000047', 'GCST006475',
        'GCST006940', 'GCST90012112', 'GCST90002398', 'GCST90002351',
        'GCST90012102', 'GCST90038600', 'GCST006943', 'GCST004624',
        'GCST90002374', 'GCST005038', 'GCST004627', 'GCST90038599',
        'GCST004617', 'GCST004600', 'GCST004606', 'GCST004623',
        'GCST004610', 'GCST011494', 'GCST006952', 'GCST90014325',
        'GCST90002407', 'GCST90002399', 'GCST90014290', 'GCST004613',
        'GCST004629', 'GCST004620', 'GCST004626', 'GCST004614',
        'GCST90002292', 'GCST004618', 'GCST004625', 'GCST004609',
        'GCST004608', 'GCST004988', 'GCST006661', 'GCST006414',
        'GCST90012106', 'GCST90012107', 'GCST90002380', 'GCST90002379',
        'GCST90012026', 'GCST011365', 'GCST90000514', 'GCST006250',
        'GCST006572', 'GCST90000048', 'GCST90012103', 'GCST90016669',
        'GCST006586', 'GCST90014023', 'GCST010681', 'GCST009722',
        'GCST004633', 'GCST004632', 'GCST006945', 'GCST006478',
        'GCST005232', 'GCST006465', 'GCST006464', 'GCST009760',
        'GCST007518', 'GCST007517', 'GCST006867', 'GCST90012114',
        'GCST002223', 'GCST002216', 'GCST90038633', 'GCST90012113',
        'GCST005413', 'GCST90014033', 'GCST004364', 'GCST90012041',
        'GCST004137', 'GCST003156', 'GCST003155', 'GCST003401',
        'GCST003372', 'GCST90038690', 'GCST90012104', 'GCST005569',
        'GCST90013410', 'GCST90014291', 'GCST003837', 'GCST90038684',
        'GCST90038683', 'GCST004631', 'GCST004634', 'GCST005061',
        'GCST90000050', 'GCST005350', 'GCST007557', 'GCST005345',
        'GCST005346', 'GCST90038656', 'GCST90038652', 'GCST007090',
        'GCST002221', 'GCST90038602', 'GCST90038679', 'GCST005527',
        'GCST90038596', 'GCST002222', 'GCST90016666', 'GCST90012002',
        'GCST90038687', 'GCST90012020', 'GCST90012790', 'GCST003496',
        'GCST90000046', 'GCST007092', 'GCST007091', 'GCST006944',
        'GCST005065', 'GCST006701', 'GCST005068', 'GCST011364',
        'GCST90012000', 'GCST005195', 'GCST009541', 'GCST90038609',
        'GCST006697', 'GCST000998', 'GCST90012794', 'GCST002783',
        'GCST90002409', 'GCST90012054', 'GCST90038689', 'GCST005536',
        'GCST90038636', 'GCST90013791', 'GCST90000529', 'GCST90000513',
        'GCST90010277', 'GCST006620', 'GCST90019411', 'GCST90019446',
        'GCST90000045', 'GCST005902', 'GCST90038681', 'GCST005581',
        'GCST003129', 'GCST006947', 'GCST004075', 'GCST90012792',
        'GCST006702', 'GCST001475', 'GCST90019404', 'GCST006099',
        'GCST004074', 'GCST009761', 'GCST006085', 'GCST004076',
        'GCST90012034', 'GCST90012014', 'GCST006268', 'GCST006948',
        'GCST90016675', 'GCST90019423', 'GCST006951', 'GCST90010170',
        'GCST90019415', 'GCST90014292', 'GCST90000290', 'GCST90000292',
        'GCST90038680', 'GCST90027161', 'GCST90010229', 'GCST90006888',
        'GCST90012025', 'GCST006098', 'GCST005067', 'GCST011073',
        'GCST90012046', 'GCST001791', 'GCST011083', 'GCST90001390',
        'GCST006906', 'GCST006908', 'GCST90012033', 'GCST90012060',
        'GCST90016667', 'GCST90011995', 'GCST010723', 'GCST90012011',
        'GCST90019403', 'GCST90019381', 'GCST90012056', 'GCST90014267',
        'GCST90012039', 'GCST006097', 'GCST005921', 'GCST90012877',
        'GCST90012878', 'GCST002245', 'GCST006946', 'GCST90020053',
        'GCST003484', 'GCST90012057', 'GCST90012007', 'GCST004460',
        'GCST90010138', 'GCST90038664', 'GCST90038661', 'GCST90012115',
        'GCST90038648', 'GCST006100', 'GCST90038615', 'GCST90038614',
        'GCST90012048', 'GCST90016673', 'GCST90012031', 'GCST90012023',
        'GCST90014289', 'GCST90014266', 'GCST006366', 'GCST004439',
        'GCST004422', 'GCST90012016', 'GCST90038628', 'GCST90038629',
        'GCST90012037', 'GCST006862', 'GCST009763', 'GCST90019433',
        'GCST90012073', 'GCST90016668', 'GCST90010296', 'GCST90012053',
        'GCST005047', 'GCST003770', 'GCST005327', 'GCST90012051',
        'GCST90012059', 'GCST003659', 'GCST90012017', 'GCST005186',
        'GCST90060189', 'GCST001212', 'GCST005314', 'GCST90019406',
        'GCST005367', 'GCST007858', 'GCST90007527', 'GCST90007526',
        'GCST90012036', 'GCST90019456', 'GCST90060270', 'GCST90060200',
        'GCST90060254', 'GCST90060178', 'GCST90060302', 'GCST90060231',
        'GCST90060298', 'GCST90038601', 'GCST90000289', 'GCST90019438',
        'GCST006942', 'GCST90019444', 'GCST90012040', 'GCST90012009',
        'GCST90012081', 'GCST90012013', 'GCST90012105', 'GCST90020091',
        'GCST90012022', 'GCST90019427', 'GCST90019482', 'GCST004734',
        'GCST004733', 'GCST90013422', 'GCST90012027', 'GCST90012050',
        'GCST90012070', 'GCST90012024', 'GCST90038619', 'GCST005058',
        'GCST005073', 'GCST90060239', 'GCST90060284', 'GCST90060308',
        'GCST90060306', 'GCST90060266', 'GCST90012080', 'GCST005647',
        'GCST90010165', 'GCST003724', 'GCST90019473', 'GCST011075',
        'GCST011078', 'GCST003375', 'GCST90012055', 'GCST006698',
        'GCST90010224', 'GCST90000515', 'GCST90012063', 'GCST90012044',
        'GCST90019384', 'GCST90010310', 'GCST90019385', 'GCST90019432',
        'GCST90010304', 'GCST90012049', 'GCST90010169', 'GCST90016670',
        'GCST004432', 'GCST90012078', 'GCST90010332', 'GCST90010195',
        'GCST004433', 'GCST90012052', 'GCST90012004', 'GCST004424',
        'GCST90012006', 'GCST90010135', 'GCST90019440', 'GCST009971',
        'GCST005923', 'GCST010005', 'GCST90007322', 'GCST90010262',
        'GCST90010150', 'GCST90016674', 'GCST90010260', 'GCST90019443',
        'GCST90010252', 'GCST90010350', 'GCST90010197', 'GCST90010334',
        'GCST90010105']

    # Args
    out_todo="liftover_commands.txt"

    todo_h=open(out_todo, "ab")

    for Study_ID in Studies_to_run:
# python Liftover_hg38_to_hg19_single_study.py 
#   --Study_ID GCSTxxxx 
#   --input_credset ~/Credible_SNP_sets/finemapping_220401.parquet/
#   --Enrichment_outdir ~/CHEERS/Results/
#   --input_peak ~/Subsampled_epimap_H3K27ac/Subsampled_Epimap_H3K27ac_counts_normToMax_quantileNorm_euclideanNorm.txt
#   --outdir ~/Credible_SNP_sets/Formatted_hg19/

        cmd=['python', "Liftover_hg38_to_hg19_single_study.py ", "--Study_ID", Study_ID, "--input_credset ~/Credible_SNP_sets/finemapping_220401.parquet/", "--Enrichment_outdir ~/CHEERS/Results/", "--input_peak ~/Subsampled_epimap_H3K27ac/Subsampled_Epimap_H3K27ac_counts_normToMax_quantileNorm_euclideanNorm.txt", "--outdir ~/Credible_SNP_sets/Formatted_hg19/"]
        cmd_str=" ".join([str(arg) for arg in cmd])
        print(cmd_str)
        todo_h.write((cmd_str + '\n').encode())
        
if __name__ == '__main__':

    main()