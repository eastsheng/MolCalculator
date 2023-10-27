# pip install pyinstaller
# conda activate py3.12
# pyinstaller -D -w ./MolCalculator.py -p ./scripts/cal_fraction.py  -p ./scripts/FindCompounds.py.py --clean  --add-data "./data; data"  --add-data "./imgs; imgs"  -i ./imgs/icons-64.png
pyinstaller MolCalculator.spec