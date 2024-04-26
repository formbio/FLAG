echo "pulling flag_augustus_singularity"
singularity pull flag_augustus.image docker://ghcr.io/formbio/flag_augustus_singularity:latest
mv flag_augustus.image containers/augustus/

echo "pulling flag_cbbasic"
singularity pull flag_cbbasic.image docker://ghcr.io/formbio/flag_cbbasic:latest
mv flag_cbbasic.image containers/cbbasic/

echo "pulling flag_combinefilter_singularity"
singularity pull flag_combinefilter.image docker://ghcr.io/formbio/flag_combinefilter_singularity:latest
mv flag_combinefilter.image containers/combinefilter/

echo "pulling flag_entap"
singularity pull flag_entap.image docker://ghcr.io/formbio/flag_entap:latest
mv flag_entap.image containers/entap/

echo "pulling flag_exonerate"
singularity pull flag_exonerate.image docker://ghcr.io/formbio/flag_exonerate:latest
mv flag_exonerate.image containers/exonerate/

# echo "pulling flag_helixer"
# singularity pull flag_helixer.image docker://ghcr.io/formbio/flag_helixer:latest
# mv flag_helixer.image containers/helixer/

echo "pulling flag_helixercpu"
singularity pull flag_helixercpu.image docker://ghcr.io/formbio/flag_helixercpu:latest
mv flag_helixercpu.image containers/helixercpu/

echo "pulling flag_liftoff"
singularity pull flag_liftoff.image docker://ghcr.io/formbio/flag_liftoff:latest
mv flag_liftoff.image containers/liftoff/

echo "pulling flag_ncbiclibraries"
singularity pull flag_ncbiclibraries.image docker://ghcr.io/formbio/flag_ncbiclibraries:latest
mv flag_ncbiclibraries.image containers/ncbiclibraries/

# echo "pulling flag_cbbasic"
# singularity pull flag_ncbitools.image docker://ghcr.io/formbio/flag_ncbitools:latest
# mv flag_ncbitools.image containers/ncbitools/

echo "pulling flag_pasa"
cp containers/augustus/flag_augustus.image containers/pasa/flag_pasa.image
# singularity pull flag_pasa.image docker://ghcr.io/formbio/flag_pasa:latest
# mv flag_pasa.image containers/pasa/

echo "pulling flag_tetools"
singularity pull flag_tetools.image docker://ghcr.io/formbio/flag_tetools:latest
mv flag_tetools.image containers/tetools/

echo "pulling flag_transdecoder"
singularity pull flag_transdecoder.image docker://ghcr.io/formbio/flag_transdecoder:latest
mv flag_transdecoder.image containers/transdecoder/

echo "pulling flag_trinity"
singularity pull flag_trinity.image docker://ghcr.io/formbio/flag_trinity:latest
mv flag_trinity.image containers/trinity/