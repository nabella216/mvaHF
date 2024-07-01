## Training

```
# ./run_TMVAClassification.sh [train] [draw curves]
./run_TMVAClassification.sh 1 1
```

## Create MVA tree for an input file
- In `run_TMVAClassification.sh`
```
inputms=(
    Input_files_MVA_will_be_applied_on.root
)
outputmva= # somewhere you can write in
```
The MVA calculation will be applied on `$inputm`
- Run
```
./run_TMVAClassification.sh 0 0 1
```
- A new tree will be added
```
   (mva_output_20240508hltZDC1n/filelist_May07_trainD0_20240508hltZDC1n_sideband_BDT_2_8_0-1-2-5-6-7-12-13_1bin.root)
   ./
    ├── (TDirectoryFile) => Dfinder
    │   ├── (TTree) => ntDkpi (556543)
    │   └── (TTree) => ntDkpi (556543)
    ├── ...
    └── (TDirectoryFile) => dataset
        ├── (TTree) => tmvainfo (1)
        └── (TTree) => mva (556543)
******************************************************************************
*Tree    :mva       :                                                        *
*Entries :   556543 : Total =        28673213 bytes  File  Size =   16351930 *
*        :          : Tree compression factor =   1.75                       *
******************************************************************************
*Br    0 :mvaDsize  : mvaDsize/I                                             *
*Entries :   556543 : Total  Size=    2233069 bytes  File Size  =     509723 *
*Baskets :       70 : Basket Size=      32000 bytes  Compression=   4.38     *
*............................................................................*
*Br    1 :BDT       : BDT[mvaDsize]/F                                        *
*Entries :   556543 : Total  Size=   19806668 bytes  File Size  =   13268696 *
*Baskets :      676 : Basket Size=      32000 bytes  Compression=   1.49     *
*............................................................................*
*Br    2 :mvapref   : mvapref[mvaDsize]/O                                    *
*Entries :   556543 : Total  Size=    6633291 bytes  File Size  =    2561718 *
*Baskets :      278 : Basket Size=      32000 bytes  Compression=   2.59     *
*............................................................................*

```

