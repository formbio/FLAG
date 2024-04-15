echo "Pulling flag_augustus"
singularity pull flag_augustus.image docker://ghcr.io/formbio/flag_augustus_singularity:latest
mv flag_augustus.image containers/augustus/

echo "Pulling flag_cbbasic"
singularity pull flag_cbbasic.image docker://ghcr.io/formbio/flag_cbbasic:latest
mv flag_cbbasic.image containers/cbbasic/

echo "Building flag_combinefilter"
cd containers/combinefilter/
cp -r ../../scripts/ .
singularity build flag_combinefilter.image Singularity.def
cd ../../

echo "Pulling flag_entap"
singularity pull flag_entap.image docker://ghcr.io/formbio/flag_entap:latest
mv flag_entap.image containers/entap/

echo "Pulling flag_exonerate"
singularity pull flag_exonerate.image docker://ghcr.io/formbio/flag_exonerate:latest
mv flag_exonerate.image containers/exonerate/

# echo "Pulling flag_augustus"
# singularity pull flag_helixer.image docker://ghcr.io/formbio/flag_helixer:latest
# mv flag_helixer.image containers/helixer/

echo "Pulling flag_helixercpu"
singularity pull flag_helixercpu.image docker://ghcr.io/formbio/flag_helixercpu:latest
mv flag_helixercpu.image containers/helixercpu/

echo "Pulling flag_liftoff"
singularity pull flag_liftoff.image docker://ghcr.io/formbio/flag_liftoff:latest
mv flag_liftoff.image containers/liftoff/

echo "Pulling flag_ncbiclibraries"
singularity pull flag_ncbiclibraries.image docker://ghcr.io/formbio/flag_ncbiclibraries:latest
mv flag_ncbiclibraries.image containers/ncbiclibraries/
# singularity pull flag_ncbitools.image docker://ghcr.io/formbio/flag_ncbitools:latest
# mv flag_ncbitools.image containers/ncbitools/

echo "Pulling flag_pasa"
singularity pull flag_pasa.image docker://ghcr.io/formbio/flag_pasa:latest
mv flag_pasa.image containers/pasa/

echo "Pulling flag_tetools"
singularity pull flag_tetools.image docker://ghcr.io/formbio/flag_tetools:latest
mv flag_tetools.image containers/tetools/

echo "Pulling flag_transdecoder"
singularity pull flag_transdecoder.image docker://ghcr.io/formbio/flag_transdecoder:latest
mv flag_transdecoder.image containers/transdecoder/

echo "Pulling flag_trinity"
singularity pull flag_trinity.image docker://ghcr.io/formbio/flag_trinity:latest
mv flag_trinity.image containers/trinity/
# singularity pull flag_augustus.image docker://ghcr.io/formbio/flag_augustus:latest
# mv flag_augustus.image containers/augustus/
