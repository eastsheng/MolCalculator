# pip install pyinstaller
pyinstaller -D -w ./MolCalculator.py --clean -i ./imgs/icons-64.png
cp -r ./imgs/ ./dist/MolCalculator/
cp -r ./scripts/ ./dist/MolCalculator/
cp -r ./temp/ ./dist/MolCalculator/
