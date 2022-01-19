python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py preproc -c /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py mapping -c /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -c -jcs 2
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py collapse -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py refine -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -refine
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py find_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py visual_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -g Zm00001d050245
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py rank_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py go -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg /data/CrazyHsu_data/Projects/isoseq/iFLAS_toolkit_20210510/iFLAS_toolkit/test_data/gene2go.txt -s sample1,sample2,sample3
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py report -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
