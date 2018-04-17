#!/bin/bash

python invariant1.py -in Hap1_1000000_iced_data.npz -out Hap1
python invariant2.py -in Hap1_1000000_iced_data.npz -out Hap1
python invariant3.py -in Hap1_1000000_iced_data.npz -out Hap1 -d 1 10
python invariant3b.py -in Hap1_1000000_iced_data.npz -out Hap1 -d 1 10
python misassembly.py -in Hap1_1000000_iced_data.npz -out misassembly -lim 750 880 -b 0.1 0.6 0.8
python misassembly.py -in Hap1_1000000_iced_data.npz -out misassembly -lim 750 880 -b 0.1 0.6 0.8 -ds 100000

