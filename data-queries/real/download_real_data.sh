#Download Real Datasets: 
echo "Downloading The Seismic Data--Start"
python  $PROJECTS_ROOT/tools/download_data/download_data.py 1jowTIVJYKqvE1qMlg1rIqNb54-P6zCvz $PROJECTS_ROOT/data/real/data_size100M_seismic_len256_znorm.bin
echo "Downloading The Seismic Data--End"

echo "Downloading The Sald Data--Start"
python  $PROJECTS_ROOT/tools/download_data/download_data.py 1bUPrKJGAVbDcWeZ5ilxnmVnnsx1t7lpG $PROJECTS_ROOT/data/real/data_size899M_sald_len128_znorm.bin
echo "Downloading The Sald Data--End"

echo "Downloading The Sift1B Data--Start"
python  $PROJECTS_ROOT/tools/download_data/download_data.py 1kWoKRMyaW2jLmOyD-xkaopAbp3GfHbJy $PROJECTS_ROOT/data/real/data_size1B_sift_len128_znorm.bin
echo "Downloading The Sift1B Data--End"

echo "Downloading The Deep1B Data--Start"
python  $PROJECTS_ROOT/tools/download_data/download_data.py 1FrSgXjg2ye4S1TxEWTqlO9uRBKxh41tw $PROJECTS_ROOT/data/real/data_size1B_deep1b_len96_znorm.bin
echo "Downloading The Deep1B Data--End"


#In case the script does not work properly, here are the shareable google drivelinks
#Seismic: https://drive.google.com/drive/u/1/folders/1jowTIVJYKqvE1qMlg1rIqNb54-P6zCvz
#Sald: https://drive.google.com/drive/u/1/folders/1bUPrKJGAVbDcWeZ5ilxnmVnnsx1t7lpG
#Deep1B: https://drive.google.com/drive/u/1/folders/1FrSgXjg2ye4S1TxEWTqlO9uRBKxh41tw
#Sift1B: https://drive.google.com/file/d/1kWoKRMyaW2jLmOyD-xkaopAbp3GfHbJy/view?usp=sharing

