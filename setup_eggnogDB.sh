# Make the eggnogDB folder
mkdir eggnogDB
cd eggnogDB
# Download the eggnog files
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz && echo Decompressing... && gunzip eggnog_proteins.dmnd.gz

wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && echo Decompressing... && gunzip eggnog.db.gz

wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && echo Decompressing... && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz

# tar the folder
cd ..
tar -zcvf eggnogDB.tar.gz eggnogDB/
