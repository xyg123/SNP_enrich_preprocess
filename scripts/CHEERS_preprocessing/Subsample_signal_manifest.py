import os
import sys
import json
import argparse
import gzip
import yaml
import pandas as pd
def main():
    # Args
    out_todo="generate_signal_commands.txt"
    Test_samples=['BSS01668', 'BSS01665', 'BSS00038', 'BSS01666', 'BSS01669', 'BSS01696', 'BSS01690', 'BSS00183', 'BSS01693', 'BSS01420', 'BSS01154', 'BSS00330', 'BSS00084', 'BSS01397', 'BSS00705', 'BSS00078', 'BSS00077', 'BSS00131', 'BSS01272', 'BSS00127', 'BSS00025', 'BSS00748', 'BSS00026', 'BSS00159', 'BSS01105', 'BSS01548', 'BSS01594', 'BSS01657', 'BSS01282', 'BSS01602', 'BSS00273', 'BSS01264', 'BSS00287', 'BSS01612', 'BSS01366', 'BSS00716', 'BSS00277', 'BSS00478', 'BSS00315', 'BSS00483', 'BSS00057', 'BSS00050', 'BSS00284', 'BSS01399', 'BSS00048', 'BSS00258', 'BSS01206', 'BSS00143', 'BSS00260', 'BSS00387', 'BSS01505', 'BSS00355', 'BSS00075', 'BSS01068', 'BSS01103', 'BSS01503', 'BSS00329', 'BSS01504', 'BSS00328', 'BSS01502', 'BSS00546', 'BSS00238', 'BSS00098', 'BSS00550', 'BSS00237', 'BSS00495', 'BSS01127', 'BSS00087', 'BSS01508', 'BSS00520', 'BSS01078', 'BSS01136', 'BSS01513', 'BSS01482', 'BSS01533', 'BSS00554', 'BSS01158', 'BSS01170', 'BSS00553', 'BSS00511', 'BSS01871', 'BSS01527', 'BSS01147', 'BSS01187', 'BSS01188', 'BSS00428', 'BSS00452', 'BSS00404', 'BSS00439', 'BSS00472', 'BSS01315', 'BSS01301', 'BSS01846', 'BSS01324', 'BSS01298', 'BSS01574', 'BSS01573', 'BSS01338', 'BSS01155', 'BSS01576', 'BSS00146', 'BSS00145', 'BSS00304', 'BSS00368', 'BSS01156', 'BSS01621', 'BSS01840', 'BSS01617', 'BSS01842', 'BSS01613', 'BSS00123', 'BSS00122', 'BSS00124', 'BSS00121', 'BSS00758', 'BSS01855', 'BSS01856', 'BSS01443', 'BSS01859', 'BSS01867', 'BSS01456', 'BSS01886', 'BSS01884', 'BSS01887', 'BSS01459', 'BSS01286', 'BSS01660', 'BSS01606', 'BSS01287', 'BSS01659', 'BSS01625', 'BSS01628', 'BSS01634', 'BSS01631', 'BSS01630', 'BSS00339', 'BSS00332', 'BSS01661', 'BSS01583', 'BSS00063', 'BSS01823', 'BSS01827', 'BSS01824', 'BSS01825', 'BSS01826', 'BSS00244', 'BSS00734', 'BSS01108', 'BSS01107', 'BSS00738']
    
    todo_h=open(out_todo, "wb")

    for sample_ID in Test_samples:
# python Generate_single_signal.py --Sample "BSS01668" --Peaks "../../tmp/Master_enhancers_d250.sorted.merged.bed" --outdir "../../tmp/H3K27ac" 
        cmd=['python', "Generate_single_signal.py", "--Sample", sample_ID, "--Peaks", "../../tmp/Master_enhancers_d250.sorted.merged.bed", "--outdir", "../../tmp/H3K27ac"]
        cmd_str=" ".join([str(arg) for arg in cmd])
        print(cmd_str)
        todo_h.write((cmd_str + '\n').encode())


if __name__ == '__main__':

    main()