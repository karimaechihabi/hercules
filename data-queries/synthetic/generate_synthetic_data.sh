GENERATOR=$PROJECTS_ROOT/tools/Cgenerator/generator/generator

#DATASETS
export GSL_RNG_SEED=1184
#sanity test
#$GENERATOR --size 1000000 --length 256  --z-normalize --filename "$DATASETS/data_size1GB_seed1184_len256_znorm.bin"

#$GENERATOR --size 2000000000 --length 256  --z-normalize --filename "$DATASETS/data_size2B_seed1184_len256_znorm.bin"
#$GENERATOR --size 200000000 --length 128     --z-normalize --filename "$DATASETS/data_size200M_seed1184_len128_znorm.bin"
#$GENERATOR --size 50000000  --length 512     --z-normalize --filename "$DATASETS/data_size50M_seed1184_len512_znorm.bin"
#$GENERATOR --size 25000000  --length 1024    --z-normalize --filename "$DATASETS/data_size25M_seed1184_len1024_znorm.bin"
#$GENERATOR --size 12500000  --length 2048    --z-normalize --filename "$DATASETS/data_size12M500K_seed16627_len2048_znorm.bin"
#$GENERATOR --size 6250000   --length 4096    --z-normalize --filename "$DATASETS/data_size6M250K_seed1184_len4096_znorm.bin"
#$GENERATOR --size 3125000   --length 8192    --z-normalize --filename "$DATASETS/data_size3M125K_seed1184_len8192_znorm.bin"
#$GENERATOR --size 1562500   --length 16384   --z-normalize --filename "$DATASETS/data_size1M562K500_seed1184_len16384_znorm.bin"

#QUERIES
#export GSL_RNG_SEED=14784
#$GENERATOR --size 1000       --length 128     --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len128_znorm.bin"
#$GENERATOR --size 1000       --length 256     --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len256_znorm.bin"
#$GENERATOR --size 1000       --length 512     --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len512_znorm.bin"
#$GENERATOR --size 1000       --length 1024    --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len1024_znorm.bin"
#$GENERATOR --size 1000       --length 2048    --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len2048_znorm.bin"
#$GENERATOR --size 1000       --length 4096    --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len4096_znorm.bin"
#$GENERATOR --size 1000       --length 8192    --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len8192_znorm.bin"
#$GENERATOR --size 1000       --length 16384   --z-normalize --filename "$DATASETS/queries_size1K_seed14784_len16384_znorm.bin"









